// std
#include <iostream>
#include <vector>
#include <fstream>
// kdtree++
#include <kdtree++/kdtree.hpp>
// teem
#include <teem/hest.h>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>
// xavier
#include <math/MLS.hpp>
#include <image/nrrd_wrapper.hpp>
#include <math/matrix.hpp>

#include <sstream>

const double particle_radius = 0.001;
const size_t max_nb_neighbors = 1000;

template<typename V>
class point {
public:
    typedef V                       vector_type;
    typedef typename V::value_type  value_type;
    
    point() : __v(), __idx(-1) {}
    point(value_type x, value_type y, value_type z, unsigned int idx = -1)
        : __v(x, y, z), __idx(idx) {}
    point(const vector_type& v, unsigned int idx = -1) : __v(v), __idx(idx) {}
    point(const point& p) : __v(p.__v), __idx(p.__idx) {}
    
    const vector_type& pos() const {
        return __v;
    }
    
    vector_type& pos() {
        return __v;
    }
    
    unsigned int index() const {
        return __idx;
    }
    
    unsigned int& index() {
        return __idx;
    }
    
    value_type distance_to(const point& p) const {
        return norm(p.__v - __v);
    }
    
    value_type operator[](size_t N) const {
        return __v[N];
    }
    
private:
    vector_type     __v;
    unsigned int    __idx;
};

template<typename V>
inline bool operator==(const point<V>& p0, const point<V>& p1)
{
    return (p0.pos() == p1.pos());
}

template<int N>
inline nvis::fixed_vector<float, N> to_vector(float* ptr)
{
    nvis::fixed_vector<float, N> r;
    for (int i = 0 ; i < N ; ++i) {
        r[i] = ptr[i];
    }
    return r;
}

inline double sane(double v)
{
    if (std::isinf(v) || std::isnan(v)) {
        return 0;
    } else {
        return v;
    }
}

char* base_in, *base_out;
int t0, dt;
float radius, dx;
void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "base input",           airTypeString,  1, 1, &base_in,             NULL,       "base input coordinates file name");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &base_out,            NULL,       "output stress field NRRD file");
    hestOptAdd(&hopt, "t0",     "initial time",         airTypeInt,     1, 1, &t0,                  NULL,       "initial time step for advection");
    hestOptAdd(&hopt, "dt",     "delta T",              airTypeInt,     1, 1, &dt,                  NULL,       "number of time steps for advection");
    hestOptAdd(&hopt, "dx",     "sampling distance",    airTypeFloat,   1, 1, &dx,                  NULL,       "distance between samples (isotropic)");
    hestOptAdd(&hopt, "r",      "support radius",       airTypeFloat,   1, 1, &radius,              "5",        "radius of weighting function support");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute per particle FTLE in granular flow",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline std::string zeros(int n)
{
    if (n < 10) {
        return std::string("000");
    } else if (n < 100) {
        return std::string("00");
    } else if (n < 1000) {
        return std::string("0");
    } else {
        return std::string("");
    }
}

inline double ftle(const Eigen::MatrixXd& jacobian)
{
    xavier::mat3 M;
    for (int i = 0 ; i < 3 ; ++i) {
        for (int j = 0 ; j < 3 ; ++j) {
            M(i, j) = jacobian(j+1, i);
        }
    }
    
    xavier::mat3 T(xavier::transpose(M));
    M *= T;
    std::vector<double> evals;
    std::vector<nvis::vec3> evecs;
    M.eigensystem(evals, evecs);
    if (evals.size()) {
        return std::max((double)0, log(*std::max_element(evals.begin(), evals.end())));
    } else {
        return 0.;
    }
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    dx *= particle_radius;
    
    std::ostringstream os;
    os << base_in << zeros(t0) << t0 << ".nrrd";
    std::string start_f = os.str();
    
    os.clear();
    os.str("");
    os << base_in << zeros(t0 + dt) << t0 + dt << ".nrrd";
    std::string fwd_f = os.str();
    
    os.clear();
    os.str("");
    os << base_in << zeros(t0 - dt) << t0 - dt << ".nrrd";
    std::string bwd_f = os.str();
    
    // import particle positions at start
    Nrrd* nin_start = nrrdNew();
    if (nrrdLoad(nin_start, start_f.c_str(), NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* coord_start = (float*)nin_start->data;
    // import particle positions after forward advection
    Nrrd* nin_fwd = nrrdNew();
    if (nrrdLoad(nin_fwd, fwd_f.c_str(), NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* coord_fwd = (float*)nin_fwd->data;
    // import particle positions after backward advection
    Nrrd* nin_bwd = nrrdNew();
    if (nrrdLoad(nin_bwd, bwd_f.c_str(), NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* coord_bwd = (float*)nin_bwd->data;
    
    // we will be working in float precision
    typedef nvis::fixed_vector<float, 3>    vector_type;
    typedef point<vector_type>              point_type;
    typedef nvis::bounding_box<vector_type> box_type;
    typedef nvis::fixed_vector<float, 3>    value_type;
    
    // construct kd-tree
    KDTree::KDTree<3, point_type> tree;
    box_type box;
    
    int nbp = nin_start->axis[1].size;
    
    // insert start positions in tree
    for (int i = 0 ; i < nbp ; ++i) {
        vector_type v(to_vector<3>(&coord_start[3*i]));
        box.add(v);
        tree.insert(point_type(v, i));
    }
    
    // now loop over a number of random positions
    vector_type min = box.min();
    vector_type max = box.max();
    vector_type diagonal = box.size();
    
    std::cout << "bounding box in input: " << min << " - " << max << std::endl;
    double cylinder_radius = max[0] + particle_radius;
    double minz = min[2] - particle_radius;
    double height = diagonal[2] + 2 * particle_radius;
    std::cerr << "radius = " << radius
              << ", minz = " << minz
              << ", height = " << height << '\n';
    box.min() -= nvis::vec3(particle_radius, particle_radius, particle_radius);
    box.max() += nvis::vec3(particle_radius, particle_radius, particle_radius);
    
    // allocate volume
    size_t nx, ny, nz;
    nx = ny = floor(2 * cylinder_radius / dx) + 1;
    nz = floor(height / dx) + 1;
    float* data = (float*)calloc(nx * ny * nz, sizeof(float));
    
    nvis::timer timer;
    double search_time = 0;
    double mls_time = 0;
    
    float distance = (float)radius * particle_radius;
    
    MLS::weighted_least_squares<value_type, nvis::vec3> WLS(3, 1, 3);
    
    for (int k = 0 ; k < nz ; ++k) {
        for (int j = 0 ; j < ny ; ++j) {
            for (int i = 0 ; i < nx ; ++i) {
            
                vector_type x = box.min() + (dx * vector_type(i, j, k));
                if (x[0]*x[0] + x[1]*x[1] > cylinder_radius*cylinder_radius) {
                    continue;
                }
                
                point_type p(x[0], x[1], x[2]);
                
                std::vector<point_type> in_cube;
                
                Eigen::MatrixXd coef_fwd(4, 3), coef_bwd(4, 3);
                std::vector<nvis::vec3> pos;
                std::vector<value_type> val_fwd, val_bwd;
                
                timer.restart();
                tree.find_within_range(p, distance, std::back_inserter(in_cube));
                search_time += timer.elapsed();
                
                if (in_cube.empty()) {
                    std::cerr << "there are no particles within distance " << radius
                              << " of " << x << std::endl;
                    continue;
                } else {
                    // std::cerr << "found " << in_cube.size() << " neighbors of " << x
                    //           << " within manhattan distance " << radius << std::endl;
                    for (int i = 0 ; i < in_cube.size() ; ++i) {
                        const point_type& close_p = in_cube[i];
                        vector_type dist = close_p.pos() - p.pos();
                        double d = nvis::norm(dist);
                        if (d > distance) {
                            continue;
                        }
                        // std::cerr << "\tparticle #" << close_p.index() << ": " << close_p.pos() << std::flush;
                        // std::cerr << " is at distance " << d << std::endl;
                        
                        pos.push_back(close_p.pos());
                        val_fwd.push_back(to_vector<3>(&coord_fwd[3*close_p.index()]));
                        val_bwd.push_back(to_vector<3>(&coord_bwd[3*close_p.index()]));
                    }
                }
                
                if (pos.size() < 4) {
                    std::cerr << "not enough particles in neighborhood (" << pos.size() << ")\n";
                    continue;
                }
                
                timer.restart();
                int deg;
                deg = WLS(coef_fwd, pos, val_fwd, x, distance);
                deg = WLS(coef_bwd, pos, val_bwd, x, distance);
                mls_time += timer.elapsed();
                
                double ftle_f = sane(ftle(coef_fwd));
                double ftle_b = sane(ftle(coef_bwd));
                double max_ftle = (ftle_f > ftle_b) ? ftle_f : -ftle_b;
                data[i + nx*(j + ny*k)] = max_ftle;
            }
        }
    }
    
    std::cerr << "Number of ftle values computed: " << 2*nx* ny* nz << '\n';
    std::cerr << "Total search time: " << search_time << " s.\n";
    std::cerr << "Average search time per sample: " << search_time / (double)nbp << " s.\n";
    std::cerr << "Total MLS time: " << mls_time << " s.\n";
    std::cerr << "Average MLS time per sample: " << mls_time / (double)nbp << " s.\n";
    
    os.clear();
    os.str("");
    os << base_out << zeros(t0) << t0 << "-ftle-dt=" << zeros(dt) << dt << ".nrrd";
    size_t dims[3] = {nx, ny, nz};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, data, nrrdTypeFloat, 3, dims);
    nrrdSave(os.str().c_str(), nout, NULL);
    
    std::cerr << "output NRRD file " << os.str() << "exported\n";
    
    return 0;
}













































































