#include <iostream>

#include <teem/hest.h>
#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/wall_timer.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in_scale;
char* name_in_val;
char* name_out;

void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me = argv[0];
    
    mop = airMopNew();
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "is",     "input scale",      airTypeString,  1,  1,  &name_in_scale,     NULL,   "input scale file name");
    hestOptAdd(&hopt, "iv",     "input value",      airTypeString,  1,  1,  &name_in_val,       NULL,   "input value file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,          NULL,   "output name");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Select value at best scale",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef xavier::nrrd_data_traits<Nrrd*>  field_type;

int main(int argc, const char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nins = xavier::nrrd_utils::readNrrd(name_in_scale);
    Nrrd* ninv = xavier::nrrd_utils::readNrrd(name_in_val);
    
    int M = nins->axis[0].size;
    int N = nins->axis[1].size;
    int nscales = nins->axis[2].size;
    int nbpix = M*N;
    
    std::vector<double> scales, values;
    xavier::nrrd_utils::to_vector(scales, nins);
    xavier::nrrd_utils::to_vector(values, ninv);
    
    double* out = (double*)calloc(M*N, sizeof(double));
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < M*N ; ++n) {
        
            double maxs = scales[n];
            int bests = 0;
            for (int s=1 ; s<nscales ; ++s) {
                double str = scales[n+s*nbpix];
                if (str>maxs) {
                    maxs = str;
                    bests = s;
                }
            }
            out[n] = values[n+bests*nbpix];
        }
    }
    
    std::cout << "done" << std::endl;
    
    std::vector<size_t> dims(2);
    dims[0] = M;
    dims[1] = N;
    std::vector<double> spc(2);
    spc[0] = nins->axis[0].spacing;
    spc[1] = nins->axis[1].spacing;
    xavier::nrrd_utils::writeNrrdFromContainers(out, name_out, /*nrrdTypeDouble,*/ dims, spc);
    
    return 0;
}


























