#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <fstream>
#include <sstream>
#include <math/fixed_vector.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* dir, *prefix, *suffix, *output;
int delta;
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
    hestOptAdd(&hopt, "f",      "folder",           airTypeString,  1, 1, &dir,     NULL,       "path of input NRRD file names");
    hestOptAdd(&hopt, "p",      "prefix",           airTypeString,  1, 1, &prefix,  NULL,       "prefix of input NRRD file names");
    hestOptAdd(&hopt, "s",      "suffix",           airTypeString,  1, 1, &suffix,  NULL,       "suffix of input NRRD file names");
    hestOptAdd(&hopt, "d",      "delta index",      airTypeInt,     0, 1, &delta,   "10",       "index step size");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &output,  NULL,       "output base name");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Determine per-pixel best iteration scale in map FTLE image",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    std::ostringstream os;
    os << dir << "/" << prefix << delta << suffix;
    
    Nrrd* nin = nrrdNew();
    nin = xavier::readNrrd(os.str());
    
    int M = nin->axis[1].size;
    int N = nin->axis[2].size;
    float hx = nin->axis[1].spacing;
    float hy = nin->axis[2].spacing;
    
    nrrdNuke(nin);
    
    int* best_scale = (int*)calloc(M * N, sizeof(int));
    float* max_str = (float*)calloc(M * N, sizeof(float));
    float* result = (float*)calloc(3 * M * N, sizeof(float));
    
    for (int it = delta ; it <= 100 ; it += delta) {
        std::ostringstream os1, os2;
        os1 << dir << '/' << prefix << it << suffix;
        os2 << dir << '/' << "lambda2-" << prefix << it << suffix;
        
        Nrrd* nin1 = nrrdNew();
        nin1 = xavier::readNrrd(os1.str());
        float* ftle = (float*)nin1->data;
        
        Nrrd* nin2 = nrrdNew();
        nin2 = xavier::readNrrd(os2.str());
        float* lambda2 = (float*)nin2->data;
        
        #pragma omp parallel for
        for (int i = 0 ; i < M*N ; ++i) {
            float str = -lambda2[i];
            if (str > max_str[i]) {
                max_str[i] = str;
                best_scale[i] = it;
                result[3*i  ] = ftle[3*i  ];
                result[3*i+2] = ftle[3*i+2];
            }
        }
    }
    
    Nrrd* nout = nrrdNew();
    os.clear();
    os.str("");
    os << output << "-best_ftle.nrrd";
    size_t dims[3] = {3, M, N};
    if (nrrdWrap_nva(nout, result, nrrdTypeFloat, 3, dims)) {
        std::cout << "ERROR while wrapping data: " << biffGetDone(NRRD)
                  << std::endl;
        if (result) {
            delete[] result;
        } else {
            nrrdNuke(nout);
        }
        exit(-1);
    }
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, airNaN(), hx, hy);
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                  << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    std::cout << "exported " << os.str() << std::endl;
    
    nout = nrrdNew();
    os.clear();
    os.str("");
    os << output << "-best_scale.nrrd";
    size_t _dims[2] = {M, N};
    if (nrrdWrap_nva(nout, best_scale, nrrdTypeFloat, 2, _dims)) {
        std::cout << "ERROR while wrapping data: " << biffGetDone(NRRD)
                  << std::endl;
        if (best_scale) {
            delete[] best_scale;
        } else {
            nrrdNuke(nout);
        }
        exit(-1);
    }
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterCell, nrrdCenterCell);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, hx, hy);
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                  << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    std::cout << "exported " << os.str() << std::endl;
    
    return 0;
}



