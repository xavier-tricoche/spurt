#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <random>
#include <math.h>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <format/DLRreader.hpp>
#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <util/timer.hpp>
#include <vtk/vtk_utils.hpp>

// parameters
size_t        nsamples = 0;
double        factor = -1;
std::string   in_name, in_name2, out_name;
nvis::bbox3   bbox;
std::vector<std::string> types;

struct Coords {
    Coords() : m_procedural(false) {}
    Coords(const Nrrd* nin) {
        initialize(nin);
    }

    Coords(VTK_SMART(vtkDataSet) dataset) {
        initialize(dataset);
    }

    void initialize(const Nrrd* nin, int first_dim=0) {
        spurt::nrrd_utils::nrrd_traits traits(nin);
        for (int d=0; d<traits.dim()-first_dim; ++d) {
            m_spacings[d] = traits.spacings()[first_dim+d];
            m_mins[d] = traits.spacings()[first_dim+d];
            m_sizes[d] = traits.sizes()[first_dim+d];
            std::cout << "spacing: " << m_spacings[d] << ", min: " << m_mins[d] << ", size: " << m_sizes[d] << '\n';
        }
        if (traits.dim() - first_dim == 2) {
            m_spacings[2] = m_mins[2] = 0;
            m_sizes[2] = 1;
        }
        m_procedural = true;
    }

    void initialize(VTK_SMART(vtkDataSet) dataset) {
        int npts = dataset->GetNumberOfPoints();
        m_pos.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            dataset->GetPoint(i, m_pos[i].begin());
        }
        m_procedural = false;
    }

    void initialize(const std::vector<nvis::vec3>& pos) {
        m_pos.resize(pos.size());
        std::copy(pos.begin(), pos.end(), m_pos.begin());
        m_procedural = false;
    }

    nvis::vec3 operator[](int i) const {
        if (!m_procedural) {
            return m_pos[i];
        }
        else {
            nvis::ivec3 id;
            id[0] = i%m_sizes[0];
            int n = i/m_sizes[0];
            if (m_spacings[2] == 0) {
                id[1] = n;
                id[2] = 0;
            }
            else {
                id[1] = n%m_sizes[1];
                id[2] = n/m_sizes[1];
            }
            return nvis::vec3(m_mins[0]+m_spacings[0]*id[0],
                              m_mins[1]+m_spacings[1]*id[1],
                              m_mins[2]+m_spacings[2]*id[2]);
        }
    }

    size_t size() const {
        if (m_procedural) {
            return m_sizes[0]*m_sizes[1]*m_sizes[2];
        }
        else {
            return m_pos.size();
        }
    }

    bool m_procedural;
    nvis::vec3 m_spacings;
    nvis::vec3 m_mins;
    nvis::ivec3 m_sizes;
    std::vector<nvis::vec3> m_pos;
};

bool has_type(const std::string& str) {
    for (int i=0; i<types.size(); ++i) {
        if (types[i] == str) return true;
    }
    return false;
}

void printUsageAndExit( const std::string& argv0, const std::string& offending="",
                        bool doExit = true )
{
    if (offending != "") {
        std::cerr << "ERROR: " << offending << std::endl;
    }
    std::cerr
        << "Usage  : " << argv0 << " [parameters] [options]\n"
        << "Synopsis: resample given dataset using smooth reconstruction"
        << "Parameters:\n"
        << " -i  | --input <string> [x2]   Input file name(s)\n"
        << " -o  | --output <string>       Output file name\n"
        << " -n  | --number <int>          Number of output samples\n"
        << " -t  | --type <string>         Value type (scalar, vector, tensor)\n"
        << "Options:\n"
        << " -b  | --bounds <float> (x6)   Bounding box of region (in world coordinates). Default: all\n"
        << std::endl;

    if (doExit) exit(1);
}

Coords all_points;

// only some of these arrays will be filled
// depending on the nature of the input data
std::vector<nvis::vec3> all_vectors;
std::vector<double>     all_scalars;
std::vector<nvis::mat3> all_tensors;

void load_VTK(const std::string& name, const std::string& me) {
    VTK_SMART(vtkDataSet) dataset = vtk_utils::readVTK(name);
    all_points.initialize(dataset);

    int npts = dataset->GetNumberOfPoints();
    VTK_SMART(vtkDataArray) scalars = dataset->GetPointData()->GetScalars();
    if (scalars != NULL && has_type("scalar")) {
        all_scalars.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_scalars[i] = *scalars->GetTuple(i);
        }
    }
    VTK_SMART(vtkDataArray) vectors = dataset->GetPointData()->GetVectors();
    if (vectors != NULL && has_type("vector")) {
        all_vectors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            vectors->GetTuple(i, all_vectors[i].begin());
        }
    }
    VTK_SMART(vtkDataArray) tensors = dataset->GetPointData()->GetTensors();
    if (tensors != NULL && has_type("tensor")) {
        double t[9];
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            tensors->GetTuple(i, t);
            int n=0;
            for (int r=0 ; r<3 ; ++r) {
                for (int c=0 ; c<3 ; ++c)
                    all_tensors[i][r][c] = t[n++];
            }
        }
    }
}

void load_NRRD(const std::string& name, const std::string& me) {
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, name.c_str(), NULL)) {
        printUsageAndExit(me, biffGetDone(NRRD));
    }
    std::vector<double> data;
    spurt::nrrd_utils::to_vector<double>(data, nin);
    // identify data type based on number of columns
    int ncol = nin->axis[0].size;
    int npts = nin->axis[1].size;
    if (nin->dim == 2 && has_type("scalar")) { // 2d scalar field
        all_points.initialize(nin, 0);
        npts = nin->axis[0].size * nin->axis[1].size;
        all_scalars.resize(npts);
        for (int i=0; i<npts; ++i) {
            all_scalars[i] = data[i];
        }
    }
    else if (nin->dim == 3 && nin->axis[0].size == 1 && has_type("scalar")) { // 2d scalar field
        all_points.initialize(nin, 1);
        npts = nin->axis[1].size * nin->axis[2].size;
        all_scalars.resize(npts);
        for (int i=0; i<npts; ++i) {
            all_scalars[i] = data[i];
        }
    }
    else if (nin->dim == 3 && nin->axis[0].size <= 3 && has_type("vector")) { // 2/3d vector field
        all_points.initialize(nin, 1);
        npts = nin->axis[1].size * nin->axis[2].size;
        all_vectors.resize(npts);
        int stride = nin->axis[0].size;
        for (int i=0 ; i<npts ; ++i) {
            all_vectors[i][0] = data[stride*i  ];
            all_vectors[i][1] = data[stride*i+1];
            if (stride == 3) all_vectors[i][2] = data[stride*i+2];
        }
    }
    else if (nin->dim == 3 && nin->axis[0].size == 4) {
        if (has_type("tensor")) { // 2d tensor field
            all_points.initialize(nin, 1);
            npts = nin->axis[1].size * nin->axis[2].size;
            all_tensors.resize(npts);
            for (int i=0; i<npts; ++i) {
                all_tensors[i][0][0] = data[4*i  ];
                all_tensors[i][0][1] = data[4*i+1];
                all_tensors[i][0][2] = 0;
                all_tensors[i][1][0] = data[4*i+2];
                all_tensors[i][1][1] = data[4*i+3];
                all_tensors[i][1][2] = 0;
                all_tensors[i][2][0] = 0;
                all_tensors[i][2][1] = 0;
                all_tensors[i][2][2] = 0;
            }
        }
        else if (has_type("scalar")) {
            npts = nin->axis[1].size;
            all_scalars.resize(npts);
            std::vector<nvis::vec3> pts(npts);
            for (int i=0 ; i<npts ; ++i) {
                for (int j=0 ; j<3 ; ++j) {
                    pts[i][j]  = data[4*i+j];
                }
                all_scalars[i] = data[4*i+3];
            }
            all_points.initialize(pts);
        }
        else {
            printUsageAndExit("Input NRRD file incompatible with specified type(s)");
        }
    }
    else if ((nin->dim == 3 || (nin->dim == 4 && nin->axis[0].size == 1)) && has_type("scalar")) {
        all_points.initialize(nin, nin->dim == 4 ? 1 : 0);
        npts = all_points.size();
        all_scalars.resize(npts);
        for (int i=0; i<npts; ++i) {
            all_scalars[i] = data[i];
        }
    }
    else if (nin->dim == 4 && nin->axis[0].size == 3 && has_type("vector")) {
        all_points.initialize(nin, 1);
        npts = all_points.size();
        all_vectors.resize(npts);
        for (int i=0; i<npts; ++i) {
            all_vectors[i][0] = data[3*i  ];
            all_vectors[i][1] = data[3*i+1];
            all_vectors[i][2] = data[3*i+2];
        }
    }
    else if (nin->dim == 4 && nin->axis[0].size == 6 && has_type("tensor")) {
        all_points.initialize(nin, 1);
        npts = all_points.size();
        all_tensors.resize(npts);
        for (int i=0; i<npts; ++i) {
            all_tensors[i][0][0] = data[6*i  ];
            all_tensors[i][0][1] = data[6*i+1];
            all_tensors[i][0][2] = data[6*i+2];
            all_tensors[i][1][0] = data[6*i+1];
            all_tensors[i][1][1] = data[6*i+3];
            all_tensors[i][1][2] = data[6*i+4];
            all_tensors[i][2][0] = data[6*i+2];
            all_tensors[i][2][1] = data[6*i+4];
            all_tensors[i][2][2] = data[6*i+5];
        }
    }
    else if (nin->dim == 4 && nin->axis[0].size == 9 && has_type("tensor")) {
        all_points.initialize(nin, 1);
        npts = all_points.size();
        all_tensors.resize(npts);
        for (int i=0; i<npts; ++i) {
            all_tensors[i][0][0] = data[6*i  ];
            all_tensors[i][0][1] = data[6*i+1];
            all_tensors[i][0][2] = data[6*i+2];
            all_tensors[i][1][0] = data[6*i+3];
            all_tensors[i][1][1] = data[6*i+4];
            all_tensors[i][1][2] = data[6*i+5];
            all_tensors[i][2][0] = data[6*i+6];
            all_tensors[i][2][1] = data[6*i+7];
            all_tensors[i][2][2] = data[6*i+8];
        }
    }
    else if (nin->dim == 2 && ncol == 6 && has_type("vector")) {
        std::vector<nvis::vec3> pts(npts);
        all_vectors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            for (int j=0 ; j<3 ; ++j) {
                pts[i][j]  = data[6*i+j];
                all_vectors[i][j] = data[6*i+3+j];
            }
        }
        all_points.initialize(pts);
    }
    else if (nin->dim == 2 && ncol == 9 && has_type("tensor")) {
        std::vector<nvis::vec3> pts(npts);
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            pts[i][0]     = data[9*i  ];
            pts[i][1]     = data[9*i+1];
            pts[i][2]     = data[9*i+2];
            all_tensors[i][0][0] = data[9*i+3];
            all_tensors[i][0][1] = data[9*i+4];
            all_tensors[i][1][0] = data[9*i+4];
            all_tensors[i][0][2] = data[9*i+5];
            all_tensors[i][2][0] = data[9*i+5];
            all_tensors[i][1][1] = data[9*i+6];
            all_tensors[i][1][2] = data[9*i+7];
            all_tensors[i][2][1] = data[9*i+7];
            all_tensors[i][2][2] = data[9*i+8];
        }
        all_points.initialize(pts);
    }
    else if (nin->dim == 2 && ncol == 12 && has_type("tensor")) {
        std::vector<nvis::vec3> pts(npts);
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            pts[i][0]     = data[12*i   ];
            pts[i][1]     = data[12*i+ 1];
            pts[i][2]     = data[12*i+ 2];
            all_tensors[i][0][0] = data[12*i+ 3];
            all_tensors[i][0][1] = data[12*i+ 4];
            all_tensors[i][0][2] = data[12*i+ 5];
            all_tensors[i][1][0] = data[12*i+ 6];
            all_tensors[i][1][1] = data[12*i+ 7];
            all_tensors[i][1][2] = data[12*i+ 8];
            all_tensors[i][2][0] = data[12*i+ 9];
            all_tensors[i][2][1] = data[12*i+10];
            all_tensors[i][2][1] = data[12*i+11];
        }
        all_points.initialize(pts);
    }
    else {
        std::cerr << "Invalid NRRD input file\n";
        exit(1);
    }
}

void load_DLR(const std::string& grid_name, const std::string data_name, const std::string& me) {
    spurt::DLRreader reader(grid_name, data_name);
    std::vector<nvis::fvec3> vertices;
    std::vector<long int> cell_indices;
    std::vector<std::pair<spurt::DLRreader::cell_type, long int> >cell_types;
    reader.read_mesh(false, vertices, cell_indices, cell_types);
    int npts = vertices.size();
    std::vector<nvis::vec3> pts(npts);
    for (int i=0 ; i<npts ; ++i) pts[i] = vertices[i];

    std::vector<double> tmp;
    if (has_type("vector")) {
        all_vectors.resize(npts);
        reader.read_data("x_velocity", tmp);
        for (int i=0 ; i<npts ; ++i) all_vectors[i][0] = tmp[i];
        reader.read_data("y_velocity", tmp);
        for (int i=0 ; i<npts ; ++i) all_vectors[i][1] = tmp[i];
        reader.read_data("z_velocity", tmp);
        for (int i=0 ; i<npts ; ++i) all_vectors[i][2] = tmp[i];
    }
    if (has_type("scalar")) {
        all_scalars.resize(npts);
        reader.read_data("pressure", tmp);
        for (int i=0 ; i<npts ; ++i) all_scalars[i] = tmp[i];
    }
}

nvis::bbox3 bounds() {
    nvis::bbox3 bb;
    for (int i=0 ; i<all_points.size() ; ++i) {
        bb.add(all_points[i]);
    }
    return bb;
}

template<typename T>
T get_val(T* array, const nvis::ivec3& id, const nvis::ivec3& size) {
    return array[id[0] + size[0]*(id[1] + size[1]*id[2])];
}

std::string base(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(0, found);
}

std::string extension(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(found+1);
}

int main(int argc, char* argv[]) {
    in_name = "none";
    in_name2 = "none";
    out_name = "none";
    bbox.min() = bbox.max() = nvis::vec3(0);

    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing input");
            }
            in_name = argv[++i];
            if (i < argc-1) {
                std::string tmp(argv[i+1]);
                if (tmp[0] != '-') {
                    in_name2 = tmp;
                    ++i;
                }
            }
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing output");
            }
            out_name = argv[++i];
        }
        else if (arg == "-n" || arg == "--number") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing number of samples");
            }
            std::string _n = argv[++i];
            if (_n[0] == 'x') {
                std::string _f = _n.substr(1);
                factor = std::stod(_f);
                std::cout << "factor=" << factor << '\n';
            }
            else if (_n[_n.size()-1] == '%') {
                std::string _f = _n.substr(0, _n.size()-1);
                int k = std::stoi(_f);
                factor = (float)k/100.;
                std::cout << "factor=" << factor << '\n';
            }
            else {
                nsamples = std::stoi(_n);
            }
        }
        else if (arg == "-t" || arg == "--type") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing type");
            }
            std::string t = argv[++i];
            spurt::lower_case(t);
            if (t != "scalar" && t != "vector" && t != "tensor") {
                printUsageAndExit(argv[0], "invalid type");
            }
            types.push_back(t);
        }
        else if (arg == "-h" || arg == "--help") {
            printUsageAndExit(argv[0]);
        }
        else {
            printUsageAndExit(argv[0], "invalid argument");
        }
    }

    if (in_name == "none" || out_name == "none") {
        printUsageAndExit(argv[0], "missing input or output file name");
    }
    else if (nsamples == 0 && factor <= 0) {
        printUsageAndExit(argv[0], "missing number of samples");
    }
    else if (types.empty()) {
        printUsageAndExit(argv[0], "missing data type(s)");
    }

    // user reader appropriate for input file type
    std::string ext = extension(in_name);
    nvis::timer _timer;
    if (ext == "vtk" || ext == "vtu" || ext == "vtp" || ext == "vtr") load_VTK(in_name, argv[0]);
    else if (ext == "nrrd") load_NRRD(in_name, argv[0]);
    else if (in_name2 != "none") {
        // DLR format with 2 files: 1 for grid and one for data
        if (ext == "grid") load_DLR(in_name, in_name2, argv[0]);
        else {
            std::string ext2 = extension(in_name2);
            if (ext2 != "grid") {
                std::ostringstream os;
                os << "unrecognized file extensions for DLR format: " << ext << " and " << ext2;
                printUsageAndExit(argv[0], os.str());
            }
            else load_DLR(in_name2, in_name, argv[0]);
        }
    }
    else {
        printUsageAndExit(argv[0], "unrecognized file type");
    }
    std::cerr << "dataset imported in " << _timer.elapsed() << " seconds\n";

    if (nvis::norm(bbox.size())) {
        nvis::bbox3 tmp = bounds();
        for (int i=0 ; i<3 ; ++i) {
            bbox.min()[i] = std::max(tmp.min()[i], bbox.min()[i]);
            bbox.max()[i] = std::min(tmp.max()[i], bbox.max()[i]);
        }
    }
    else bbox = bounds();
    std::cout << "bounding box = " << bbox << std::endl;

    size_t nrhs = 0;
    if (all_scalars.size()) nrhs++;
    if (all_vectors.size()) nrhs += 3;
    if (all_tensors.size()) nrhs += 9;

    if (factor > 0) {
        nsamples = all_points.size()*factor;
    }
    std::cout << nsamples << " points to select\n";

    _timer.restart();
    std::vector<unsigned int> indices(all_points.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    nvis::bbox3 out_box;

    ext = extension(out_name);
    if (ext == "nrrd") {
        int stride = nrhs + 3;
        nvis::fvec3 *data = (nvis::fvec3 *)calloc(stride*nsamples, sizeof(nvis::fvec3));
        for (int i=0; i<nsamples ; ++i) {
            nvis::vec3 p = all_points[indices[i]];
            out_box.add(p);
            data[stride*i  ] = p[0];
            data[stride*i+1] = p[1];
            data[stride*i+2] = p[2];
            int k=3;
            if (all_scalars.size()) {
                data[stride*i+(k++)] = all_scalars[indices[i]];
            }
            if (all_vectors.size()) {
                data[stride*i+(k++)] = all_vectors[indices[i]][0];
                data[stride*i+(k++)] = all_vectors[indices[i]][1];
                data[stride*i+(k++)] = all_vectors[indices[i]][2];
            }
            if (all_tensors.size()) {
                data[stride*i+(k++)] = all_tensors[indices[i]][0][0];
                data[stride*i+(k++)] = all_tensors[indices[i]][0][1];
                data[stride*i+(k++)] = all_tensors[indices[i]][0][2];
                data[stride*i+(k++)] = all_tensors[indices[i]][1][0];
                data[stride*i+(k++)] = all_tensors[indices[i]][1][1];
                data[stride*i+(k++)] = all_tensors[indices[i]][1][2];
                data[stride*i+(k++)] = all_tensors[indices[i]][2][0];
                data[stride*i+(k++)] = all_tensors[indices[i]][2][1];
                data[stride*i+(k++)] = all_tensors[indices[i]][2][2];
            }
        }

        size_t sizes[2] = {nrhs+3, nsamples};
        spurt::nrrd_utils::writeNrrd(data, out_name, nrrdTypeDouble, 2, sizes, false);
    }
    else if (ext == "txt") {
        std::ofstream ofs;
        ofs.open(out_name);
        ofs << nsamples << ' ' << nrhs + 3 << '\n';
        for (int i=0 ; i<nsamples ; ++i) {
            nvis::vec3 p = all_points[indices[i]];
            out_box.add(p);
            ofs << p[0] << ' ' << p[1] << ' ' << p[2];
            if (all_scalars.size()) {
                ofs << ' ' << all_scalars[indices[i]];
            }
            if (all_vectors.size()) {
                ofs << ' ' << all_vectors[indices[i]][0];
                ofs << ' ' << all_vectors[indices[i]][1];
                ofs << ' ' << all_vectors[indices[i]][2];
            }
            if (all_tensors.size()) {
                ofs << ' ' << all_tensors[indices[i]][0][0];
                ofs << ' ' << all_tensors[indices[i]][0][1];
                ofs << ' ' << all_tensors[indices[i]][0][2];
                ofs << ' ' << all_tensors[indices[i]][1][0];
                ofs << ' ' << all_tensors[indices[i]][1][1];
                ofs << ' ' << all_tensors[indices[i]][1][2];
                ofs << ' ' << all_tensors[indices[i]][2][0];
                ofs << ' ' << all_tensors[indices[i]][2][1];
                ofs << ' ' << all_tensors[indices[i]][2][2];
            }
            ofs << '\n';
        }
        ofs.close();
    }
    else if (ext == "vtk") {
        // do something
        std::cout << "we should be doing something but we are not\n";
    }

    std::cout << "bounding box on output: " << out_box << '\n';
    std::cerr << out_name << " exported\n";
    return 0;
}
