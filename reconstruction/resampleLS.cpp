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
#include <Eigen/Core>
#include <Eigen/SVD>
#include <util/timer.hpp>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkPointSet.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <sfcnn.hpp>
#include <map>

// parameters
nvis::ivec3        resolution;
std::string        in_name, in_name2, out_name;
nvis::bbox3        bbox;

void printUsageAndExit( const std::string& argv0, const std::string& offending="", 
                        bool doExit = true )
{
    if (offending != "") {
        std::cerr << "ERROR: " << offending << std::endl;
    }
      std::cerr 
    << "Usage  : " << argv0 << " [parameters] [options]\n"
    << "Synopsis: resample given dataset through global least squares\n"
    << "Parameters:\n"
    << "    -i  | --input <string> [x2]     Input file name(s)\n"
    << "    -o  | --output <string>         Output file name\n"
    << "    -r  | --resolution <int> (x3)   Resampling resolution in X, Y, and Z\n"
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
    xavier::nrrd_utils::to_vector<double>(data, nin);
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
    xavier::DLRreader reader(grid_name, data_name);
    std::vector<nvis::fvec3> vertices;
    std::vector<long int> cell_indices;
    std::vector<std::pair<xavier::DLRreader::cell_type, long int> >cell_types;
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

std::string base(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(0, found);
}

std::string extension(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(found+1);
}

struct entry {
    int i, j;
    float value;
};
void interpolate(std::vector<entry>& A, int n, const nvis::vec3& index) {
    const nvis::ivec3& r = resolution;
    static const int offset[] = {
        0,                // (0,0,0)
        1,                // (1,0,0)
        r[0]+1,           // (1,1,0)
        r[0],             // (0,1,0)
        r[0]*r[1],        // (0,0,1)
        r[0]*r[1]+1,      // (1,0,1)
        r[0]*r[1]+r[0]+1, // (1,1,1)
        r[0]*r[1]+r[0]    // (0,1,1)
    };
    
    // vertex local coordinates
    nvis::ivec3 id(floor(index[0]), floor(index[1]), floor(index[2]));
    nvis::vec3 pos = index - nvis::vec3(id);
    const double& u = pos[0];
    const double& v = pos[1];
    const double& w = pos[1];
    double weights[] = {    
        (1-u)*(1-v)*(1-w), 
        u*(1-v)*(1-w), 
        u*v*(1-w), 
        (1-u)*v*(1-w),
        (1-u)*(1-v)*w, 
        u*(1-v)*w, 
        u*v*w, 
        (1-u)*v*w
    };
    
    int base_id = id[0] + r[0]*(id[1] + r[1]*id[2]);
    for (int i=0 ; i<8 ; ++i) {
        A.push_back(entry());
        entry& e = A.back();
        e.i = n;
        e.j = base_id + offset[i];
        e.value = weights[i];
    }
}

int main(int argc, char* argv[]) {
    in_name = "none";
    in_name2 = "none";
    out_name = "none";
    resolution = nvis::ivec3(0);
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
        else {
            printUsageAndExit(argv[0], "invalid argument");
        }
    }
    
    if (in_name == "none" || out_name == "none") {
        printUsageAndExit(argv[0], "missing input or output file name");
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
    nvis::vec3 diameter = bbox.size();
    bbox.min() -= 0.001*diameter;
    bbox.max() += 0.001*diameter;
    std::cout << "bounding box = " << bbox << std::endl;
    
    nvis::vec3 spacing = bbox.size() / nvis::vec3(resolution - nvis::ivec3(1,1,1));
    std::cout << "spacing = " << spacing << std::endl;
    
    size_t nrhs = 0;
    if (all_scalars.size()) nrhs++;
    if (all_vectors.size()) nrhs += 3;
    if (all_tensors.size()) nrhs += 9;
    
    size_t nrows = all_points.size();
    size_t ncols = resolution[0]*resolution[1]*resolution[2];
    
    std::vector<float> rhs_mat(nrows*nrhs);
    
    size_t nb_threads = 1;
    
#if _OPENMP
    nb_threads = omp_get_max_threads();
#endif    
    std::cout << nb_threads << " threads available\n";
    std::vector<entry> ls_mat[nb_threads];
    
#pragma omp parallel
    {
#pragma omp for schedule (dynamic, 1) 
    for (size_t n=0 ; n<nrows ; ++n) {
        size_t thread_id = 0;
#if _OPENMP
        thread_id = omp_get_thread_num();
#endif
        
        const nvis::vec3& x = all_points[n];
        nvis::vec3 y = x - bbox.min();
        nvis::vec3 u = y / spacing;
        interpolate(ls_mat[thread_id], n, u);
        size_t offset = n*nrhs;
        if (all_scalars.size()) rhs_mat[offset++] = all_scalars[n];
        if (all_vectors.size()) {
            rhs_mat[offset++] = all_vectors[n][0];
            rhs_mat[offset++] = all_vectors[n][1];
            rhs_mat[offset++] = all_vectors[n][2];
        }
        if (all_tensors.size()) {
            for (int r=0 ; r<3 ; ++r) 
                for (int c=0 ; c<3 ; ++c)
                    rhs_mat[offset++] = all_tensors[n](r,c);
        }
        
        if (!thread_id) {
            std::cout << "\rProgress: " << std::setw(7) << std::setfill(' ') << n << " rows (" 
            << std::setw(3) << std::setfill(' ') << 100*n/nrows << "%) in "
            << _timer.elapsed() << " seconds (" << (double)n/_timer.elapsed() << " Hz)   ";
        }
    }
    }
    std::cout << std::endl;
    
    std::cout << "\nComputing reindexing... " << std::flush;
    std::map<int, int> old2new;
    size_t ndof = 0;
    for (size_t k=0 ; k<nb_threads ; ++k) {
        for (size_t i=0 ; i<ls_mat[k].size() ; ++i) {
            size_t n = ls_mat[k][i].j;
            old2new[n] = ndof++;
        }
    }
    std::cout << "done.\n";
    std::cout << "there are " << ndof << " (" << 100*ndof/ncols << "%) degrees of freedom in this problem\n";
    
    std::cout << "exporting reindexing array... " << std::flush;
    std::string ind_name;
    {
        std::ostringstream os;
        os << out_name << "_index.txt";
        ind_name = os.str();
    }
    std::fstream idf(ind_name.c_str(), std::ios::out);
    for (std::map<int, int>::const_iterator it=old2new.begin() ; it!=old2new.end() ; ++it) {
        idf << it->first << " " << it->second << '\n';
    }
    idf.close();
    
    std::cout << "exporting LS matrix..." << std::flush;
    std::string ls_name;
    {
        std::ostringstream os;
        os << out_name << "_LS.txt";
        ls_name = os.str();
    }
    std::fstream lsf(ls_name.c_str(), std::ios::out);
    for (size_t k=0 ; k<nb_threads ; ++k) {
        for (size_t i=0 ; i<ls_mat[k].size() ; ++i) {
            lsf << ls_mat[k][i].i+1 << " \t" << old2new[ls_mat[k][i].j]+1 << " \t" << ls_mat[k][i].value << '\n';
        }
    }
    lsf << nrows << " \t" << ndof << " 0\n";
    lsf.close();
    std::cout << " done\n";
    std::cout << "exporting RHS matrix..." << std::flush;
    std::string rhs_name;
    {
        std::ostringstream os;
        os << out_name << "_RHS.txt";
        rhs_name = os.str();
    }    
    std::fstream rhsf(rhs_name.c_str(), std::ios::out);
    for (size_t i=0 ; i<rhs_mat.size()/nrhs ; ++i) {
        size_t offset = i*nrhs;
        for (size_t j=0 ; j<nrhs ; ++j) {
            rhsf << rhs_mat[offset++];
            if (j<nrhs-1) rhsf << " \t";
        }
        rhsf << '\n';
    } 
    rhsf.close();
    std::cout << " done.\n";
    
    return 0;
}
