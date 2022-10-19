#include <iostream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <sstream>
#include <kdtree/kdtree.hpp>

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
int res[2], n_indices;
double dx, t, pour, relaxation, _gamma, ymin, ymax;

template<typename T, size_t N>
struct distance_traits< nvis::fixed_vector<T, N> > {
    typedef double value_type;
    
    template<unsigned int D>
    static value_type dist1( const nvis::fixed_vector<T, N>& p0,
                             const nvis::fixed_vector<T, N>& p1 ) {
        return p0[D] - p1[D];
    }
    
    static value_type dist( const nvis::fixed_vector<T, N>& p0,
                            const nvis::fixed_vector<T, N>& p1 ) {
        return nvis::norm(p0-p1);
    }
};

typedef kdtree<nvis::fvec3, int, 3> kdtree_type;

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
    hestOptAdd(&hopt, "i",      "input",                airTypeString,  1, 1, &name_in,         NULL,       "input NRRD file (3D)");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &name_out,        NULL,       "output NRRD file (2D)");
    hestOptAdd(&hopt, "r",      "resolution",           airTypeInt,     2, 2, &res,             NULL,       "output image resolution");
    hestOptAdd(&hopt, "dx",     "precision",            airTypeDouble,  0, 1, &dx,              "0.01",     "sampling step size along each ray");
    hestOptAdd(&hopt, "t",      "time step",            airTypeDouble,  0, 1, &t,               "0",        "simulation time step");
    hestOptAdd(&hopt, "p",      "pour",                 airTypeDouble,  0, 1, &pour,            "0",        "pouring time");
    hestOptAdd(&hopt, "rel",    "relaxation",           airTypeDouble,  0, 1, &relaxation,      "0",        "relaxation time");
    hestOptAdd(&hopt, "g",      "gamma",                airTypeDouble,  0, 1, &_gamma,          "0",        "tap\'s gamma coefficient");
    hestOptAdd(&hopt, "ymax",   "max height",           airTypeDouble,  0, 1, &ymax,            "0.65",     "max height in particle column");
    hestOptAdd(&hopt, "n",      "# indices",            airTypeInt,     0, 1, &n_indices,       "3456",     "total number of particle indices");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute solid fraction of particle assembly and project onto X=const plane",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

double read(kdtree_type& tree, const std::string& name, std::vector<double>& density)
{
    Nrrd* nin = xavier::readNrrd(name);
    std::vector<double> data;
    xavier::to_vector<double>(data, nin);
    int N = nin->axis[1].size;
    std::cerr << "there are " << N << " points\n";
    density.resize(n_indices);
    
    const double K = 4./3.*M_PI*(0.01*0.01*0.01)/0.740480489;
    
    double miny = std::numeric_limits<double>::max();
    for (int i=0 ; i<N ; ++i) {
        nvis::fvec3 p;
        p[0] = data[5*i+1];
        p[1] = data[5*i+2];
        p[2] = data[5*i+3];
        int id = data[5*i];
        tree.add(p, id);
        density[id] = K/data[5*i+4];
        miny = std::min(miny, data[5*i+2]);
    }
    tree.sort();
    
    std::cerr << "Data read, tree built\n";
    return miny;
}

double trace_one(const nvis::vec2& yz, kdtree_type& tree, const std::vector<double>& density)
{
    nvis::fvec3 p(-0.01, yz[0], yz[1]);
    std::vector<kdtree_type::const_iterator> n_kd(10);
    
    if (yz[0] < ymin) {
        return 0;
    }
    
    double d = 0;
    for (; p[0]<0.25; p[0]+=dx) {
        tree.find_n_nearest(p, 1, n_kd.begin());
        d += dx*density[n_kd.front()->second];
    }
    return d/0.26;
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    double floorh=0;
    ymin = 0;
    if (t > pour) {
        double bump_time = 0.5/7.5;
        double tap_time = relaxation + bump_time;
        double loc_t = fmod(t-pour, tap_time);
        
        std::cerr << "t=" << t << ", pour = " << pour << ", relaxation=" << relaxation << ", bump_time=" << bump_time
                  << ", tap_time=" << tap_time << ", loc_t=" << loc_t << std::endl;
                  
        if (loc_t < bump_time) {
            double omega = 2.*M_PI*7.5;
            double amplitude = 9.81*_gamma/(omega*omega);
            floorh = amplitude*sin(omega*loc_t);
            
            std::cerr << "omega=" << omega << ", amplitude=" << amplitude << ", flh=" << ymin << std::endl;
            std::cerr << "omega*loc_t/(2*pi) = " << omega* loc_t/2./M_PI << std::endl;
        }
    } else {
        std::cerr << "t=" << t << " < relaxation=" << relaxation << std::endl;
    }
    std::cerr << "floor height = " << floorh << std::endl;
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    std::vector<double> density;
    kdtree_type tree;
    ymin = read(tree, name_in, density);
    
    std::cerr << "column height = " << ymin << std::endl;
    ymin = floorh;
    
    nvis::bbox3 bounds(nvis::vec3(-0.01,0,-0.01), nvis::vec3(0.25,ymax,0.25));
    double dy = bounds.size()[1] / res[1];
    double dz = bounds.size()[2] / res[0];
    float* img = (float*)calloc(res[0]*res[1], sizeof(float));
    
    #pragma omp parallel for
    for (int i=0 ; i<res[0]*res[1] ; ++i) {
        int k = i%res[0];
        int j = i/res[0];
        nvis::vec2 yz(bounds.min()[1] + (float)j*dy, bounds.min()[2] + (float)k*dz);
        img[i] = trace_one(yz, tree, density);
    }
    
    std::vector<size_t> size(2);
    size[0] = res[0];
    size[1] = res[1];
    std::vector<double> spacing(2);
    spacing[0] = dz;
    spacing[1] = dy;
    xavier::writeNrrdFromContainers(img, name_out, /*nrrdTypeFloat, */size, spacing);
    
    exit(0);
}
