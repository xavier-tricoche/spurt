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
#include <data/locator.hpp>
#include <math/types.hpp>
#include <format/dlr_reader.hpp>
#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <math/mls.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <misc/progress.hpp>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkPointSet.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>


#include <vtk/vtk_utils.hpp>
#include <misc/progress.hpp>

using namespace spurt;

// parameters
svec3  resolution;
std::string  in_name, in_name2, out_name;
int          degree;
double       radius;
bbox3        bbox;

void printUsageAndExit( const std::string& argv0, 
                        const std::string& offending="",
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
        << " -r  | --resolution <int> (x3) Resampling resolution in X, Y, and Z\n"
        << " -rd | --radius <float>        Radius of MLS fit (in index space)\n"
        << "Options:\n"
        << " -d  | --degree <int>          Polynomial degree of reconstruction. Default: 0\n"
        << " -b  | --bounds <float> (x6)   Bounding box of region (in world coordinates). Default: all\n"
        << std::endl;

    if (doExit) exit(1);
}

typedef spurt::point_locator<vec3, int>   locator_type;

std::vector<vec3> all_points;
// only some of these arrays will be filled
// depending on the nature of the input data
std::vector<vec3>   all_vectors;
std::vector<double> all_scalars;
std::vector<mat3>   all_tensors;

void load_VTK(const std::string& name, const std::string& me) {
    VTK_SMART(vtkDataSet) dataset = vtk_utils::readVTK(name);
    int npts = dataset->GetNumberOfPoints();
    all_points.resize(npts);
    for (int i=0 ; i<npts ; ++i) {
        dataset->GetPoint(i, &all_points[i][0]);
    }
    VTK_SMART(vtkDataArray) scalars = dataset->GetPointData()->GetScalars();
    if (scalars != NULL) {
        all_scalars.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_scalars[i] = *scalars->GetTuple(i);
        }
    }
    VTK_SMART(vtkDataArray) vectors = dataset->GetPointData()->GetVectors();
    if (vectors != NULL) {
        all_vectors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            vectors->GetTuple(i, &all_vectors[i][0]);
        }
    }
    VTK_SMART(vtkDataArray) tensors = dataset->GetPointData()->GetTensors();
    if (tensors != NULL) {
        double t[9];
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            tensors->GetTuple(i, t);
            int n=0;
            for (int r=0 ; r<3 ; ++r) {
                for (int c=0 ; c<3 ; ++c)
                    all_tensors[i](r,c) = t[n++];
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
            all_tensors[i](0,0) = data[9*i+3];
            all_tensors[i](0,1) = data[9*i+4];
            all_tensors[i](1,0) = data[9*i+4];
            all_tensors[i](0,2) = data[9*i+5];
            all_tensors[i](2,0) = data[9*i+5];
            all_tensors[i](1,1) = data[9*i+6];
            all_tensors[i](1,2) = data[9*i+7];
            all_tensors[i](2,1) = data[9*i+7];
            all_tensors[i](2,2) = data[9*i+8];
        }
    }
    else if (ncol == 12) {
        all_points.resize(npts);
        all_tensors.resize(npts);
        for (int i=0 ; i<npts ; ++i) {
            all_points[i][0]     = data[12*i   ];
            all_points[i][1]     = data[12*i+ 1];
            all_points[i][2]     = data[12*i+ 2];
            all_tensors[i](0,0) = data[12*i+ 3];
            all_tensors[i](0,1) = data[12*i+ 4];
            all_tensors[i](0,2) = data[12*i+ 5];
            all_tensors[i](1,0) = data[12*i+ 6];
            all_tensors[i](1,1) = data[12*i+ 7];
            all_tensors[i](1,2) = data[12*i+ 8];
            all_tensors[i](2,0) = data[12*i+ 9];
            all_tensors[i](2,1) = data[12*i+10];
            all_tensors[i](2,1) = data[12*i+11];
        }
    }
    else {
        std::cerr << "Invalid NRRD input file\n";
        exit(1);
    }
}

void load_DLR(const std::string& grid_name, const std::string data_name, 
              const std::string& me) {
    spurt::dlr_reader reader(grid_name, data_name);
    std::vector<fvec3> vertices;
    std::vector<long int> cell_indices;
    std::vector<std::pair<spurt::dlr_reader::cell_type, long int> >cell_types;
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

bbox3 bounds() {
    bbox3 bb;
    for (int i=0 ; i<all_points.size() ; ++i) {
        bb.add(all_points[i]);
    }
    return bb;
}

template<int N>
Eigen::VectorXd do_fit(const vec3& x0,
    const std::vector<int>& neighs, double radius, int nrhs, int ndof) {

    typedef small_vector<double, N> val_type;
    typedef Eigen::VectorXd res_type;
    typedef spurt::MLS::weighted_least_squares<val_type, vec3> fit_type;

    std::vector<vec3> points(neighs.size());
    std::vector<val_type> values(neighs.size());
    fit_type fit(3, degree, nrhs);

    for (int i=0 ; i<neighs.size() ; ++i) {
        int id = neighs[i];
        points[i] = all_points[id];
        int c=0;
        if (all_scalars.size()) values[i][c++] = all_scalars[id];
        if (all_vectors.size()) {
            for (int k=0 ; k<3 ; ++k) {
                values[i][c++] = all_vectors[id][k];
            }
        }
        if (all_tensors.size()) {
            for (int row=0 ; row<3 ; ++row) {
                for (int col=0 ; col<3 ; ++col) {
                    values[i][c++] = all_tensors[id](row,col);
                }
            }
        }
    }

    Eigen::MatrixXd m;
    int precision = fit(m, points, values, x0, radius);
    res_type res(ndof*nrhs);
    int i=0;
    for (int r=0 ; r<m.rows() ; ++r) {
        for (int c=0 ; c<m.cols() ; ++c) {
            res(i++) = m(r,c);
        }
    }
    for (; i<ndof*nrhs ; ++i) res(i) = 0;

    return res;
}

template<typename Locator>
void which_neighbors(std::vector<int>& ns, const vec3& x0, double r, 
                     Locator& pl) {
    typedef Locator                                 locator_type;
    typedef typename locator_type::point_type       point_type;
    typedef typename locator_type::const_iterator   const_iterator;

    std::list<point_type> neighbors;
    ns.clear();
    pl.find_within_range(neighbors, x0, r);
    typename std::list<point_type>::const_iterator it;
    for (it=neighbors.begin(); it!=neighbors.end() ; ++it) {
        if (norm(it->position()-x0) < r) {
            ns.push_back(it->data());
        }
    }
}

template<typename Locator>
int how_many(const vec3& x0, double r, Locator& pl) {
    std::vector<int> ns;
    which_neighbors<Locator>(ns, x0, r, pl);
    return ns.size();
}

template<typename Locator>
double what_radius(const vec3& x0, int& N, double minrad, double maxrad, 
                   double minstep, Locator& pl) {
    // adaptively determine needed sampling radius to solve WLS problem
    double lo_r, hi_r;
    lo_r = minrad;
    hi_r = 0;

    // use radius doubling to obtain an upper bound on the radius
    int count;
    for (double r=minrad ; r<maxrad ; r=std::min(maxrad, 2*r)) {
        count = how_many<Locator>(x0, r, pl);
        if (count >= N) {
            hi_r = r;
            break;
        }
        else {
            lo_r = r;
        }
    }
    if (hi_r == 0) {
        N = count;
        return maxrad;
    }

    // use binary search to arrive at min needed radius
    while (hi_r - lo_r > minstep) {
        double r = 0.5*(lo_r + hi_r);
        int count = how_many<Locator>(x0, r, pl);
        if (count >= N) hi_r = r;
        else lo_r = r;
    }

    N = count;
    return hi_r;
}

template<typename Locator>
double what_radius(const vec3& x0, int N, Locator& nnl) {
    std::vector<unsigned long> ids;
    std::vector<double> dist;
    auto nn = nnl.find_n_nearest(x0, N);
    double maxd = 0;
    std::for_each(nn.begin(), nn.end(), [&](auto p) {
        if (p.second > maxd) maxd = p.second;
    });
    return maxd;
}

template<typename T>
T get_val(T* array, const ivec3& id, const svec3& size) {
    return array[id[0] + size[0]*(id[1] + size[1]*id[2])];
}

void low_pass(float* radii, const svec3& res, int niter, float alpha=0.5) {
    int npts = res[0]*res[1]*res[2];
    float* before = (float *)calloc(npts, sizeof(float));
    float* after = (float *)calloc(npts, sizeof(float));
    for (int i=0 ; i<npts ; ++i) before[i] = radii[i];

    int r01 = res[0]*res[1];
    const int& r0 = res[0];

    // apply umbrella operator
    for (int n=0 ; n<niter ; ++n) {
#pragma openmp parallel for
        for (int l=0 ; l<npts ; ++l) {
            int k = l/r01;
            int m = l%r01;
            int j = m/r0;
            int i = m%r0;
            int c = 0;
            float sum = 0;
            if (i>0) {
                ivec3 id(i-1, j, k);
                sum += get_val<float>(before, id, res);
                ++c;
            }
            if (i<res[0]-1) {
                ivec3 id(i+1, j, k);
                sum += get_val<float>(before, id, res);
                ++c;
            }
            if (j>0) {
                ivec3 id(i, j-1, k);
                sum += get_val<float>(before, id, res);
                ++c;
            }
            if (j<res[1]-1) {
                ivec3 id(i, j+1, k);
                sum += get_val<float>(before, id, res);
                ++c;
            }
            if (k>0) {
                ivec3 id(i, j, k-1);
                sum += get_val<float>(before, id, res);
                ++c;
            }
            if (k<res[2]-1) {
                ivec3 id(i, j, k+1);
                sum += get_val<float>(before, id, res);
                ++c;
            }
            sum /= (float)c;
            after[l] = alpha*sum + (1-alpha)*before[l];
        }
    }

    // determine needed offset to restore minimum radius at each vertex
    float offset = 0;
    for (int i=0 ; i<npts ; ++i) {
        float dr = radii[i] - after[i];
        if (dr > offset) offset = dr;
    }
    for (int i=0 ; i<npts ; ++i) {
        radii[i] = after[i] + offset;
    }

    free(before);
    free(after);
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
    radius = 0;
    degree = 0;
    resolution = spurt::svec3(0);
    bbox.min() = bbox.max() = vec3(0);

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
        else if (arg == "-d" || arg == "--degree") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing degree");
            }
            degree = atoi(argv[++i]);
        }
        else if (arg == "-rd" || arg == "--radius") {
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
    ProgressDisplay timer(false);
    if (ext == "vtk" || ext == "vtu") load_VTK(in_name, argv[0]);
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
    std::cerr << "dataset imported in " << timer.wall_time() << " seconds\n";

    if (norm(bbox.size())) {
        bbox3 tmp = bounds();
        for (int i=0 ; i<3 ; ++i) {
            bbox.min()[i] = std::max(tmp.min()[i], bbox.min()[i]);
            bbox.max()[i] = std::min(tmp.max()[i], bbox.max()[i]);
        }
    }
    else bbox = bounds();
    std::cout << "bounding box = " << bbox << std::endl;

    vec3 spacing = bbox.size() / resolution;
    std::cout << "spacing = " << spacing << std::endl;

    typedef spurt::point_locator<vec3, int>      locator_type;
    typedef locator_type::point_type             point_type;
    typedef locator_type::region_type            region_type;

    locator_type locator;
    timer.start();
    for (int i=0 ; i<all_points.size() ; ++i) {
        locator.insert(point_type(all_points[i], i));
    }
    std::cout << "point locator created in " << timer.wall_time() << " seconds\n";

    double rad = radius * *std::max_element(spacing.begin(), spacing.end());
    double maxrad = 0.5*norm(bbox.size());
    double step = 0.1*rad;
    size_t nrhs = 0;
    if (all_scalars.size()) nrhs++;
    if (all_vectors.size()) nrhs += 3;
    if (all_tensors.size()) nrhs += 9;
    int ndof = spurt::MLS::dof(3, degree);
    std::cout << "there are " << ndof << " DOF per scalar attribute for polynomial order " << degree << '\n';

    size_t npts = spurt::product(resolution);
    std::cout << npts << " sample points to evaluate\n";
    size_t r01 = resolution[0]*resolution[1];
    const size_t& r0 = resolution[0];

    std::cout << ndof*nrhs << " scalars will be estimated at each sample point\n";
    float* out_data = (float *)calloc(npts*ndof*nrhs, sizeof(float));

    timer.start();
    int minn = std::numeric_limits<int>::max();
    int maxn = 0;
    double mint = std::numeric_limits<double>::max();
    double maxt = 0;
    std::cout << "\nDetermining sampling scale at each vertex...\n";
    float* radii = (float *)calloc(npts, sizeof(float));

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
        vec3 x = bbox.min() + vec3(i,j,k) * spacing;
        int N = ndof;
        ProgressDisplay rad_t(false);
        rad_t.start();
        // radii[n] = what_radius(x, N, rad, maxrad, step, locator);
        radii[n] = std::max(what_radius(x, N, locator), rad);
        double dt = rad_t.wall_time();

        if (!thread_id) {
            if (N < minn) minn = N;
            if (N > maxn) maxn = N;
            if (dt < mint) mint = dt;
            if (dt > maxt) maxt = dt;
            std::cout << "\rProgress: " << std::setw(7) << std::setfill(' ') << n << " vertices ("
            << std::setw(3) << std::setfill(' ') << 100*n/npts << "%) in "
            << timer.wall_time() << " seconds (" << (double)n/timer.wall_time() << " Hz), "
            << " #neighbors: " << minn << " to " << maxn << ", query time: "
            << mint << " to " << maxt <<  "   ";
        }
    }
    float minr = *std::min_element(radii, &radii[npts]);
    float maxr = *std::max_element(radii, &radii[npts]);
    {
        // NRRD file storage
        size_t size[] = { resolution[0], resolution[1], resolution[2] };
        double spc[] = { spacing[0], spacing[1], spacing[2] };
        Nrrd *nout = nrrdNew();
        if (nrrdWrap_nva(nout, radii, nrrdTypeFloat, 3, size)) {
            std::cerr << argv[0] << ':' << biffGetDone(NRRD) << '\n';
            exit(1);
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        std::string name = base(out_name) + "-scale.nrrd";
        if (nrrdSave(name.c_str(), nout, NULL)) {
            std::cerr << argv[0] << ':' << biffGetDone(NRRD) << '\n';
            exit(1);
        }
        std::cout << "Exported " << name << std::endl;
    }

    timer.start();
    low_pass(radii, resolution, 10);
    std::cout << "low pass filtering of scale information took " << timer.wall_time() << " seconds\n";
    {
        // NRRD file storage
        size_t size[] = { resolution[0], resolution[1], resolution[2] };
        double spc[] = { spacing[0], spacing[1], spacing[2] };
        Nrrd *nout = nrrdNew();
        if (nrrdWrap_nva(nout, radii, nrrdTypeFloat, 3, size)) {
            std::cerr << argv[0] << ':' << biffGetDone(NRRD) << '\n';
            exit(1);
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        std::string name = base(out_name) + "-smoothscale.nrrd";
        if (nrrdSave(name.c_str(), nout, NULL)) {
            std::cerr << argv[0] << ':' << biffGetDone(NRRD) << '\n';
            exit(1);
        }
        std::cout << "Exported " << name << std::endl;
    }

    timer.start();
    std::cout << "\nComputing WLS solution at each vertex...\n";
    typedef std::list<point_type> list_type;
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
        vec3 x = bbox.min() + vec3(i,j,k) * spacing;

        std::vector<int> close;
        which_neighbors(close, x, radii[n], locator);

        if (!thread_id) {
            std::cout << "\rProgress: sampled " << std::setw(7) << std::setfill(' ') << n << " vertices ("
            << std::setw(3) << std::setfill(' ') << 100*n/npts << "%) in "
            << timer.wall_time() << " seconds (" << (double)n/timer.wall_time() << " Hz)  ";
        }

        // antialiasing
        if (close.size() > 100*ndof) {
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(close.begin(), close.end(), g);
            close.resize(100*ndof);
        }

        Eigen::VectorXd res;
        switch (nrhs) {
        case 1:
            res = do_fit<1>(x, close, rad, nrhs, ndof);
            break;
        case 3:
            res = do_fit<3>(x, close, rad, nrhs, ndof);
            break;
        case 4:
            res = do_fit<4>(x, close, rad, nrhs, ndof);
            break;
        case 9:
            res = do_fit<9>(x, close, rad, nrhs, ndof);
            break;
        case 10:
            res = do_fit<10>(x, close, rad, nrhs, ndof);
            break;
        case 12:
            res = do_fit<12>(x, close, rad, nrhs, ndof);
            break;
        case 13:
            res = do_fit<13>(x, close, rad, nrhs, ndof);
            break;
        default:
            std::cerr << "invalid dimension for right hand side: " << nrhs << std::endl;
            exit(1);
        }

        for (int l=0 ; l<ndof*nrhs ; ++l) {
            out_data[ndof*nrhs*n+l] = res(l);
        }
    }

    double elapsed = timer.wall_time();
    std::cout
         << "Statistics:\n"
        << "Algorithm time : " << elapsed << " seconds for " << npts << " points\n"
        << "Performance: " << 1000.*elapsed/(double)npts << " ms per sample\n";

    // NRRD file storage
    size_t size[] = { ndof*nrhs, resolution[0], resolution[1], resolution[2] };
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
