#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <format/DLRreader.hpp>
#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <math/MLS.hpp>
#include <math/RBFbasis.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <util/timer.hpp>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkPointSet.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <sfcnn.hpp>
#include <math/RBF.hpp>

// parameters
nvis::fixed_vector<size_t, 3>        resolution;
std::string        in_name, in_name2, out_name;
int             degree;
double            radius;
nvis::bbox3        bbox;

void printUsageAndExit( const std::string& argv0, const std::string& offending="",
                        bool doExit = true )
{
    if (offending != "") {
        std::cerr << "ERROR: " << offending << std::endl;
    }
      std::cerr
    << "Usage  : " << argv0 << " [parameters] [options]\n"
    << "Synopsis: resample given dataset using smooth RBF reconstruction\n"
    << "Parameters:\n"
    << "    -i  | --input <string> [x2]     Input file name(s)\n"
    << "    -o  | --output <string>         Output file name\n"
    << "    -r  | --resolution <int> (x3)   Resampling resolution in X, Y, and Z\n"
    << "    -x  | --radius <float>          Multiples of min RBF radius to be used\n"
    << "Options:\n"
    << "    -b  | --bounds <float> (x6)     Bounding box of region (in world coordinates). Default: all\n"
    << std::endl;

    if (doExit) exit(1);
}


std::vector<nvis::vec3> all_points;
// only some of these arrays will be filled
// depending on the nature of the input data
std::vector<nvis::vec3> all_vectors;
std::vector<double>     all_scalars;
std::vector<nvis::mat3> all_tensors;

void load_VTK(const std::string& name, const std::string& me) {
    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName(name.c_str());
    reader->Update();
    vtkDataSet* dataset = reader->GetOutput();
    int npts = dataset->GetNumberOfPoints();
    all_points.resize(npts);
    for (int i=0 ; i<npts ; ++i) {
        dataset->GetPoint(i, all_points[i].begin());
    }
    vtkDataArray* scalars = dataset->GetPointData()->GetScalars();
    if (scalars != NULL) {
        all_scalars.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_scalars[i] = *scalars->GetTuple(i);
        }
        scalars->Delete();
    }
    vtkDataArray* vectors = dataset->GetPointData()->GetVectors();
    if (vectors != NULL) {
        all_vectors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            vectors->GetTuple(i, all_vectors[i].begin());
        }
        vectors->Delete();
    }
    vtkDataArray *tensors = dataset->GetPointData()->GetTensors();
    if (tensors != NULL) {
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
        tensors->Delete();
    }
    dataset->Delete();
    reader->Delete();
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
    if (ncol == 4) {
        all_points.resize(npts);
        all_scalars.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_points[i][0] = data[4*i  ];
            all_points[i][1] = data[4*i+1];
            all_points[i][2] = data[4*i+2];
            all_scalars[i]   = data[4*i+3];
        }
    }
    else if (ncol == 6) {
        all_points.resize(npts);
        all_vectors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            for (int j=0 ; j<3 ; ++j) {
                all_points[i][j]  = data[6*i+j];
                all_vectors[i][j] = data[6*i+3+j];
            }
        }
    }
    else if (ncol == 9) {
        all_points.resize(npts);
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_points[i][0]     = data[9*i  ];
            all_points[i][1]     = data[9*i+1];
            all_points[i][2]     = data[9*i+2];
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
    }
    else if (ncol == 12) {
        all_points.resize(npts);
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_points[i][0]     = data[12*i   ];
            all_points[i][1]     = data[12*i+ 1];
            all_points[i][2]     = data[12*i+ 2];
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
    all_points.resize(npts);
    for (int i=0 ; i<npts ; ++i) all_points[i] = vertices[i];
    all_vectors.resize(npts);
    std::vector<double> tmp;
    reader.read_data("x_velocity", tmp);
    for (int i=0 ; i<npts ; ++i) all_vectors[i][0] = tmp[i];
    reader.read_data("y_velocity", tmp);
    for (int i=0 ; i<npts ; ++i) all_vectors[i][1] = tmp[i];
    reader.read_data("z_velocity", tmp);
    for (int i=0 ; i<npts ; ++i) all_vectors[i][2] = tmp[i];
    all_scalars.resize(npts);
    reader.read_data("pressure", tmp);
    for (int i=0 ; i<npts ; ++i) all_scalars[i] = tmp[i];
}

nvis::bbox3 bounds() {
    nvis::bbox3 bb;
    for (int i=0 ; i<all_points.size() ; ++i) {
        bb.add(all_points[i]);
    }
    return bb;
}

template<typename Locator>
double what_radius(const nvis::vec3& x0, int N, Locator& nnl) {
    std::vector<unsigned long> ids;
    std::vector<double> dist;
    nnl.ksearch(x0, N+1, ids, dist);
    return *std::max_element(dist.begin(), dist.end());
}

std::string base(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(0, found);
}

std::string extension(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(found+1);
}

struct myRBF {
    double _R;
    myRBF(double radius) : _R(radius) {}
    double radius() const { return _R; }
    double operator()(double r) const { return spurt::RBF::wendland(r, _R); }
};

int main(int argc, char* argv[]) {
    in_name = "none";
    in_name2 = "none";
    out_name = "none";
    radius = 0;
    degree = 0;
    resolution = nvis::fixed_vector<size_t, 3>(0);
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
        else if (arg == "-h" || arg == "--help") {
            printUsageAndExit(argv[0]);
        }
        else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-3) {
                printUsageAndExit(argv[0], "missing resolution information");
            }
            resolution[0] = atoi(argv[++i]);
            resolution[1] = atoi(argv[++i]);
            resolution[2] = atoi(argv[++i]);
        }
        else if (arg == "-x" || arg == "--radius") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing radius");
            }
            radius = atof(argv[++i]);
        }
        else {
            printUsageAndExit(argv[0], "invalid argument");
        }
    }

    if (in_name == "none" || out_name == "none") {
        printUsageAndExit(argv[0], "missing input or output file name");
    }
    else if (radius <= 0) {
        printUsageAndExit(argv[0], "missing / invalid radius information");
    }
    else if (*std::min_element(resolution.begin(), resolution.end()) <= 0) {
        printUsageAndExit(argv[0], "missing / invalid resolution information");
    }

    // user reader appropriate for input file type
    std::string ext = extension(in_name);
    nvis::timer _timer;
    if (ext == "vtk") load_VTK(in_name, argv[0]);
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

    nvis::vec3 spacing = bbox.size() / nvis::vec3(resolution);
    std::cout << "spacing = " << spacing << std::endl;

    typedef sfcnn<nvis::fvec3, 3, float>            NNlocator_type;

    size_t nrhs = 0;
    if (all_scalars.size()) nrhs++;
    if (all_vectors.size()) nrhs += 3;
    if (all_tensors.size()) nrhs += 9;

    int npts = resolution[0]*resolution[1]*resolution[2];
    std::cout << npts << " sample points to evaluate\n";
    int r01 = resolution[0]*resolution[1];
    const int& r0 = resolution[0];

    _timer.restart();
    nvis::fvec3 *__all_points = (nvis::fvec3 *)calloc(npts, sizeof(nvis::fvec3));
    for (int i=0 ; i<npts ; ++i) __all_points[i] = all_points[i];
    NNlocator_type nnl(__all_points, npts);
    std::cout << "nearest neighbor search data structure created in " << _timer.elapsed() << " seconds\n";

    std::cout << nrhs << " scalars will be estimated at each sample point\n";
    float* out_data = (float *)calloc(npts*nrhs, sizeof(float));

#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif

    _timer.restart();
    std::cout << "\nDetermining needed radius at each vertex...\n";
    std::vector<float> radii(npts, 0);

    const int nneeded = spurt::MLS::dof(3, 2);

#pragma omp parallel
    {
#pragma omp for schedule (dynamic, 1)
        for (int n=0 ; n<npts ; ++n) {
            int thread_id = 0;
#if _OPENMP
        thread_id = omp_get_thread_num();
#endif

            int k = n/r01;
            int m = n%r01;
            int j = m/r0;
            int i = m%r0;
            nvis::vec3 x = bbox.min() + nvis::vec3(i,j,k)*spacing;
            nvis::timer rad_t;
            radii[n] = what_radius(x, nneeded, nnl);
            double dt = rad_t.elapsed();
    
            if (!thread_id) {
                std::cout << "\rProgress: " << std::setw(7) << std::setfill(' ') << n << " vertices ("
                << std::setw(3) << std::setfill(' ') << 100*n/npts << "%) in "
                << _timer.elapsed() << " seconds (" << (double)n/_timer.elapsed() << " Hz)   ";
            }
        }
    }

    std::sort(radii.begin(), radii.end());
    std::cout << "Scale statistics: \n";
    for (int i=0 ; i<100 ; i+=5) {
        std::cout << std::setw(3) << std::setfill(' ') << i << "%: " << radii[i*npts/100];
    }
    std::cout << "100%: " << radii.back() << "\n";

    float R = radius*radii.back();

    std::vector<Eigen::VectorXd> all_values(npts);
    for (int i=0 ; i<npts ; ++i) {
        all_values[i] = Eigen::VectorXd(nrhs);
        int c=0;
        if (all_scalars.size()) all_values[i][c++] = all_scalars[i];
        if (all_vectors.size())
            for (int j=0 ; j<3 ; ++j) all_values[i][c++] = all_vectors[i][j];
        if (all_tensors.size())
            for (int row=0 ; row<3 ; ++row)
                for (int col=0 ; col<3 ; ++col) all_values[i][c++] = all_tensors[i](row,col);
    }
    myRBF phi(R);
    std::cout << "Solving RBF system...\n";
    _timer.restart();

    typedef Eigen::VectorXd     data_type;
    spurt::RBF::CompactSupportRBFInterpolator<data_type, double, 3, myRBF>
        interpolator(all_points, all_values, phi);
    std::cout << "RBF solution computed in " << _timer.elapsed() << " seconds\n";

    std::cout << "\nComputing RBF interpolation at each vertex...\n";
#pragma openmp parallel for
    for (int n=0 ; n<npts ; ++n) {
        int thread_id = 0;
#if _OPENMP
        thread_id = omp_get_thread_num();
#endif

        int k = n/r01;
        int m = n%r01;
        int j = m/r0;
        int i = m%r0;
        nvis::vec3 x = bbox.min() + nvis::vec3(i,j,k)*spacing;
        if (!thread_id) {
            std::cout << "\rProgress: sampled " << std::setw(7) << std::setfill(' ') << n << " vertices ("
            << std::setw(3) << std::setfill(' ') << 100*n/npts << "%) in "
            << _timer.elapsed() << " seconds (" << (double)n/_timer.elapsed() << " Hz)  ";
        }

        data_type res = interpolator(x);

        for (int l=0 ; l<nrhs ; ++l) {
            out_data[nrhs*n+l] = res(l);
        }
    }

    double elapsed = _timer.elapsed();
    std::cout
         << "Statistics:\n"
        << "Algorithm time : " << elapsed << " seconds for " << npts << " points\n"
        << "Performance: " << 1000.*elapsed/(double)npts << " ms per sample\n";

    // NRRD file storage
    size_t size[] = { nrhs, resolution[0], resolution[1], resolution[2] };
    double spc[] = { airNaN(), spacing[0], spacing[1], spacing[2] };
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_data, nrrdTypeFloat, 4, size)) {
        std::cerr << argv[0] << ':' << biffGetDone(NRRD) << '\n';
        exit(1);
    }
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
    if (nrrdSave(out_name.c_str(), nout, NULL)) {
        std::cerr << argv[0] << ':' << biffGetDone(NRRD) << '\n';
        exit(1);
    }
    std::cerr << out_name << " exported\n";
    return 0;
}
