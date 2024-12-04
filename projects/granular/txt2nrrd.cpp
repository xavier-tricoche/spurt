#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <teem/nrrd.h>

char* in_name, *out_base;
size_t nparticles, first, period;
void init(int argc, char* argv[])
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
    hestOptAdd(&hopt, "i",      "input file",       airTypeString,  1, 1, &in_name,     NULL,   "input file name (txt)");
    hestOptAdd(&hopt, "o",      "output base",      airTypeString,  1, 1, &out_base,    NULL,   "output base name");
    hestOptAdd(&hopt, "n",      "# particles",      airTypeInt,     0, 1, &nparticles,  "3456", "number of particles per time step");
    hestOptAdd(&hopt, "f",      "first time step",  airTypeInt,     0, 1, &first,       "0",    "first time step to consider");
    hestOptAdd(&hopt, "p",      "step period",      airTypeInt,     0, 1, &period,      "1",    "interval between exported steps");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert NJIT simulation data to NRRD format",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

bool read_line_coord(float& t, float& x, float& y, float& z, std::fstream& file)
{
    int i;
    file >> i >> x >> y >> z >> t;
    return !file.eof();
}
int main(int argc, char* argv[])
{
    init(argc, argv);
    
    std::fstream input(in_name, std::ios::in);
    if (!input) {
        std::cerr << argv[0] << ": unable to open " << in_name;
        exit(-1);
    }
    
    size_t size[2] = {3, nparticles};
    std::ostringstream os;
    
    int counter = 0;
    float* data;
    for (; !input.eof() ; ++counter) {
    
        os.clear();
        os.str("");
        os << out_base << "_" << std::setw(6) << std::setfill('0') << counter << ".nrrd";
        
        float x, y, z, t;
        int i;
        if (counter >= first && !((counter-first)%period)) {
            data = (float*)calloc(3*nparticles, sizeof(float));
        }
        
        for (i=0 ; i<nparticles && !input.eof() ; ++i) {
            if (!read_line_coord(t, x, y, z, input)) {
                break;
            }
            if (counter >= first && !((counter-first)%period)) {
                data[3*i  ] = x;
                data[3*i+1] = y;
                data[3*i+2] = z;
            }
        }
        if (first <= counter && !((counter-first)%period) && i == nparticles) {
            Nrrd* nout = nrrdNew();
            if (nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, size) ||
                    nrrdSave(os.str().c_str(), nout, NULL)) {
                std::cerr << biffGetDone(NRRD) << std::endl;
                exit(-1);
            }
            std::cerr << "exported:\t" << os.str() << ", t=" << t << " s." << std::endl;
            nrrdNuke(nout);
        } else if (i < nparticles) {
            std::cerr << "time step prematurily interrupted. exit" << std::endl;
            exit(-1);
        }
    }
    
    return 0;
    
}


















