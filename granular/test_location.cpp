#include <iostream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <sfcnn.hpp>
#include <sstream>
#include <kdtree/kdtree.hpp>
#include <util/timer.hpp>

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
int nsamples;

typedef sfcnn<nvis::fvec3, 3, float> NNtree_type;

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
    hestOptAdd(&hopt, "n",      "# samples",            airTypeInt,     1, 1, &nsamples,        NULL,       "number of samples");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "compare NN results between STANN and Christoph's Kdtree",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    int nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cout << nthreads << " threads available\n";
    
    nvis::bounding_box<nvis::fvec3> bounds(nvis::fvec3(0,0,0), nvis::fvec3(0.24,0.53,0.24));
    
    Nrrd* nin = xavier::readNrrd(name_in);
    std::vector<double> data;
    xavier::to_vector<double>(data, nin);
    int N = nin->axis[1].size;
    std::cerr << "there are " << N << " points\n";
    nvis::fvec3* all_points = new nvis::fvec3[N];
    
    kdtree_type tree2;
    for (int i=0 ; i<N ; ++i) {
        nvis::fvec3& p = all_points[i];
        p[0] = data[3*i  ];
        p[1] = data[3*i+1];
        p[2] = data[3*i+2];
        tree2.add(p, i);
    }
    NNtree_type tree1(all_points, N);
    tree2.sort();
    
    srand48(time(0));
    
    std::vector<float> timer_nn(nthreads, 0);
    std::vector<float> timer_kd(nthreads, 0);
    std::vector<int> wrong(nthreads, 0);
    
    #pragma omp parallel for
    for (int i=0 ; i<nsamples ; ++i) {
        int th = 0;
#if _OPENMP
        th = omp_get_thread_num();
#endif
        nvis::fvec3 query = bounds.min() + nvis::fvec3(drand48(), drand48(), drand48())*bounds.size();
        
        std::vector<long unsigned int> n_nn;
        std::vector<kdtree_type::const_iterator> n_kd(10);
        
        nvis::timer t;
        tree1.ksearch(query, 1, n_nn);
        timer_nn[th] += t.elapsed();
        
        t.restart();
        tree2.find_n_nearest(query, 1, n_kd.begin());
        timer_kd[th] += t.elapsed();
        
        if (n_kd.front()->second != n_nn.front()) {
            ++wrong[th];
        }
    }
    
    float time1=0, time2=0;
    int errors = 0;
    for (int i=0 ; i<nthreads ; ++i) {
        time1 += timer_nn[i];
        time2 += timer_kd[i];
        errors += wrong[i];
    }
    
    std::cout << "STANN took " << time1 << " (" << time1/(float)nsamples << " Hz)\n"
              << "Kd-tree took " << time2 << " (" << time2/(float)nsamples << " Hz)\n"
              << "there were " << errors << " errors (" << errors*100/nsamples << "%)\n";
              
    exit(0);
}
