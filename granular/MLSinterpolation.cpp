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

const double particle_radius = 0.05;
const size_t max_nb_neighbors = 100;

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

char* coord_f, *stress_f, *output_f;
int n;
float radius, h;
void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "c",      "coordinates",          airTypeString,  1, 1, &coord_f,             NULL,       "input coordinates file name");
    hestOptAdd(&hopt, "s",      "stress",               airTypeString,  1, 1, &stress_f,            NULL,       "input stress value file");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &output_f,            NULL,       "output stress field NRRD file");
    hestOptAdd(&hopt, "n",      "# samples",            airTypeInt,     0, 1, &n,                   "100000",   "number of random samples");
    hestOptAdd(&hopt, "r",      "support radius",       airTypeFloat,   0, 1, &radius,              "0.2",      "radius of weighting function support");
    hestOptAdd(&hopt, "h",      "sampling step",        airTypeFloat,   0, 1, &h,                   "0",        "sampling step on raster grid");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute smooth MLS interpolation of principal stress tensor over cylindrical domain",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, const char* argv[])
{
    initialize(argc, argv);
    
    // import particle positions
    Nrrd* nin1 = nrrdNew();
    if (nrrdLoad(nin1, coord_f, NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    // import particle stress tensor
    Nrrd* nin2 = nrrdNew();
    if (nrrdLoad(nin2, stress_f, NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    
    // we will be working in float precision
    typedef nvis::fixed_vector<float, 3>    vector_type;
    typedef point<vector_type>              point_type;
    typedef nvis::bounding_box<vector_type> box_type;
    typedef nvis::fixed_vector<float, 6>    value_type;
    
    float* coords = (float*)nin1->data;
    float* stress = (float*)nin2->data;
    
    // number of valid particles in file
    int nbp = nin2->axis[1].size; // only valid positions have associated stress tensor
    
    // construct kd-tree
    KDTree::KDTree<3, point_type> tree;
    box_type box;
    
    // insert positions in tree
    for (int i = 0 ; i < nbp ; ++i) {
        vector_type v(to_vector<3>(&coords[3*i]));
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
    double height = diagonal[2] + 0.1;
    std::cerr << "radius = " << radius
              << ", minz = " << minz
              << ", height = " << height << '\n';
    box.min() -= nvis::vec3(particle_radius, particle_radius, particle_radius);
    box.max() += nvis::vec3(particle_radius, particle_radius, particle_radius);
    
    // srand48(time(0));
    
    if (h <= 0 && n > 0) {
        h = pow(4.*cylinder_radius * cylinder_radius * height / n, 0.33333333);
    }
    unsigned int nx = (unsigned int)floor(2.*cylinder_radius / h);
    unsigned int nz = (unsigned int)floor(height / h);
    n = nx * nx * nz;
    std::cerr << "sampling step = " << h << ", number of particles = " << n << '\n';
    
    // allocate volume
    float* data = (float*)calloc(n * 7, sizeof(float));
    
    MLS::weighted_least_squares<value_type, nvis::vec3> WLS(3, 0, 6);
    
    nvis::timer timer;
    double search_time = 0;
    double mls_time = 0;
    for (int k = 0 ; k < nz ; ++k) {
        for (int j = 0 ; j < nx ; ++j) {
            for (int i = 0 ; i < nx ; ++i) {
            
                int id = i + nx * (j + nx * k);
                data[7*id] = 1;
                
                nvis::vec3 x = nvis::vec3((float)i * h, (float)j * h, (float)k * h);
                x += box.min();
                
                if (x[0]*x[0] + x[1]*x[1] > cylinder_radius*cylinder_radius) {
                    continue;
                }
                
                point_type p(x[0], x[1], x[2]);
                
                std::vector<point_type> in_cube;
                
                Eigen::MatrixXd coef(1,6); // first-order fit of symmetric matrix
                std::vector<nvis::vec3> pos;
                std::vector<value_type> val;
                
                timer.restart();
                tree.find_within_range(p, radius, std::back_inserter(in_cube));
                search_time += timer.elapsed();
                
                if (in_cube.empty()) {
                    std::cerr << "there are no particles within distance " << radius
                              << " of " << x << std::endl;
                    continue;
                } else {
                    for (int i = 0 ; i < in_cube.size() ; ++i) {
                        const point_type& close_p = in_cube[i];
                        vector_type dist = close_p.pos() - p.pos();
                        double d = nvis::norm(dist);
                        if (d > radius) {
                            continue;
                        }
                        
                        pos.push_back(close_p.pos());
                        val.push_back(to_vector<6>(&stress[close_p.index()*6]));
                    }
                }
                
                if (pos.size() < 3) {
                    std::cerr << "not enough particles in neighborhood (" << pos.size() << ")\n";
                    continue;
                }
                
                timer.restart();
                int deg = WLS(coef, pos, val, x, radius);
                mls_time += timer.elapsed();
                
                for (int j = 0 ; j < 6 ; ++j) {
                    data[7*id+j] = coef(0,j);
                }
            }
        }
    }
    
    std::cerr << "Number of tensor values computed: " << n << '\n';
    std::cerr << "Total search time: " << search_time << " s.\n";
    std::cerr << "Average search time per sample: " << search_time / (double)n << " s.\n";
    std::cerr << "Total MLS time: " << mls_time << " s.\n";
    std::cerr << "Average MLS time per sample: " << mls_time / (double)n << " s.\n";
    
    std::vector<double> spacing(4);
    spacing[0] = airNaN();
    spacing[1] = spacing[2] = spacing[3] = h;
    std::vector<size_t> size(4);
    size[0] = 7;
    size[1] = size[2] = nx;
    size[3] = nz;
    xavier::writeNrrdFromContainers(data, output_f, /*nrrdTypeFloat,*/ size, spacing);
    
    std::cerr << "output NRRD file exported\n";
    
    return 0;
}




































































