#include <teem/hest_helper.hpp>
#include <image/nrrd_wrapper.hpp>
#include <iostream>
#include <vector>


char* in, *out;
int prec;
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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &in,      NULL,       "input file name (NRRD)");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &out,     NULL,       "output file name (TIFF)");
    hestOptAdd(&hopt, "p",      "precision",        airTypeInt,     0, 1, &prec,    "4",        "encoding precision (in bytes)");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Convert 3-channel FTLE NRRD file to color map in TIFF format",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin = nrrdNew();
    nin = spurt::readNrrd(in);
    
    // verify data type
    if (nin->dim != 3 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> __array;
    spurt::to_vector(__array, nin);
    int M = nin->axis[1].size;
    int N = nin->axis[2].size;
    
    nvis::bbox2 value_range;
    for (int i = 0 ; i < M*N ; ++i) {
        nvis::vec2 val(__array[3*i], __array[3*i+2]);
        value_range.add(val);
    }
    
    std::cerr << "Value range = " << value_range  << '\n';
    
    float* colors = (float*)calloc(3 * M * N, sizeof(float));
    
    for (int i = 0 ; i < M*N ; ++i) {
        nvis::vec2 val(__array[3*i], __array[3*i+2]);
        val -= value_range.min();
        val /= value_range.size();
        double u = 1 - val[0];
        double v = 1 - val[1];
        nvis::fvec3 c = u * (1 - v) * nvis::fvec3(0, 0, 1) + v * (1 - u) * nvis::fvec3(1, 0, 0) + u * v * (1, 1, 1);
        colors[3*i  ] = c[0];
        colors[3*i+1] = c[1];
        colors[3*i+2] = c[2];
    }
    
    Nrrd* nout = nrrdNew();
    size_t dims[] = {3, M, N};
    if (nrrdWrap_nva(nout, colors, nrrdTypeFloat, 3, dims)) {
        std::cout << "ERROR while wrapping data: " << biffGetDone(NRRD)
                  << std::endl;
        if (colors) {
            delete[] colors;
        } else {
            nrrdNuke(nout);
        }
        exit(-1);
    }
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, airNaN(), nin->axis[1].spacing, nin->axis[2].spacing);
    if (nrrdSave(out, nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                  << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    std::cout << "exported " << out << std::endl;
    
    return 0;
}


