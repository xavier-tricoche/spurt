#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <sfcnn.hpp>
#include <sstream>
#include <map>
#include <iomanip>

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in, *name_out;
float height;
float epsilon;
float t, relaxation, _gamma;

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
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &name_out,        NULL,       "output TXT file");
    hestOptAdd(&hopt, "e",      "epsilon",              airTypeFloat,   0, 1, &epsilon,         "0",        "contact tolerance");
    hestOptAdd(&hopt, "t",      "time step",            airTypeFloat,   0, 1, &t,               "0",        "simulation time step");
    hestOptAdd(&hopt, "rel",    "relaxation",           airTypeFloat,   0, 1, &relaxation,      "0",        "relaxation time");
    hestOptAdd(&hopt, "g",      "gamma",                airTypeFloat,   0, 1, &_gamma,          "0",        "tap\'s gamma coefficient");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute number of contacts in particle column with periodic lateral boundaries",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef sfcnn<nvis::fvec3, 3, float> tree_type;

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    size_t nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    
    height = 0;
    if (t > relaxation) {
        float bump_time = 0.5/7.5;
        float tap_time = relaxation + bump_time;
        float loc_t = fmod(t-relaxation, tap_time);
        if (loc_t < bump_time) {
            float omega = 2.*M_PI*7.5;
            float amplitude = 9.81*_gamma/(omega*omega);
            height = amplitude*sin(omega*loc_t);
        }
    }
    std::cerr << "floor height = " << height << std::endl;
    
    Nrrd* nin = spurt::readNrrd(name_in);
    std::vector<float> data;
    spurt::to_vector<float>(data, nin);
    int N = nin->axis[1].size;
    std::vector<nvis::fvec3> points(N);
    
#pragma openmp parallel for
    for (int i=0 ; i<N ; ++i) {
        nvis::fvec3& p = points[i];
        p[0] = data[3*i  ];
        p[1] = data[3*i+1];
        p[2] = data[3*i+2];
    }
    
    const float min_c = 0;
    const float max_c = 0.24;
    const float span = max_c - min_c;
    const float radius = 0.01;
    
    typedef std::pair<nvis::fvec3, int> duplicate_type;
    std::vector<std::vector<duplicate_type> > duplicates(nthreads);
    
#pragma openmp parallel for
    for (int i=0 ; i<N ; ++i) {
        const nvis::fvec3& p = points[i];
        
        int th_id = 0;
#if _OPENMP
        th_id = omp_get_thread_num();
#endif
        
        if (p[0] < min_c+radius) {
            duplicates[th_id].push_back(duplicate_type(p + nvis::fvec3(span, 0, 0), i));
        } else if (p[0] > max_c - radius) {
            duplicates[th_id].push_back(duplicate_type(p - nvis::fvec3(span, 0, 0), i));
        }
        if (p[2] < min_c+radius) {
            duplicates[th_id].push_back(duplicate_type(p + nvis::fvec3(0, 0, span), i));
        } else if (p[2] > max_c - radius) {
            duplicates[th_id].push_back(duplicate_type(p - nvis::fvec3(0, 0, span), i));
        }
    }
    
    int nduplicates = 0;
    for (int i=0 ; i<nthreads ; ++i) {
        nduplicates += duplicates[i].size();
    }
    
    std::cerr << "there are " << nduplicates << " duplicates\n";
    nvis::fvec3* all_points = new nvis::fvec3[N+nduplicates];
    for (int i=0 ; i<N ; ++i) {
        all_points[i] = points[i];
    }
    int count = N;
    
    std::map<int, int> duplicate_to_originals;
    
    for (int n=0 ; n<nthreads ; ++n) {
        for (int i=0 ; i<duplicates[n].size() ; ++i) {
            all_points[count] = duplicates[n][i].first;
            duplicate_to_originals[count++] = duplicates[n][i].second;
        }
    }
    
    tree_type tree(all_points, N+nduplicates);
    
    int nppc = 0;
    int npfc = 0;
#pragma openmp parallel for
    for (int i=0 ; i<N ; ++i) {
        int th_id = 0;
#if _OPENMP
        th_id = omp_get_thread_num();
#endif
        
        std::vector<unsigned long> answers;
        std::vector<double> distances;
        tree.ksearch(points[i], 12, answers, distances);
        
        if (points[i][1] < height + radius + epsilon) {
            ++npfc;
        }
        
        for (int j=0 ; j<answers.size() ; ++j) {
            if (distances[j] > 2*radius + epsilon) {
                break;
            }
            unsigned long id = answers[j];
            if (id >= N) {
                id = duplicate_to_originals[id];
            }
            if (id > i) {
                ++nppc;
            }
        }
    }
    
    std::cerr << (float)nppc / (float)N << " particle contacts on average\n";
    
    std::fstream output(name_out, std::ios::out | std::ios::app);
    output << t << " " << N << " " << nppc << " " << npfc << '\n';
    output.close();
    
    exit(0);
}
