#include <iostream>
#include <iomanip>

#include <vector>
#include <complex>
#include <sstream>
#include <math.h>

#include <boost/format.hpp>
#include <boost/limits.hpp>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>

// spurt
#include <math/math.hpp>

// christoph
#include <tokamak/poincare_map.hpp>
#include <tokamak/tokamak_nimrod_parametric.hpp>

#include "definitions.hpp"

#include <teem/hest_helper.hpp>


char* file, *outs, *ts;
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
    hestOptAdd(&hopt, "f",      "file",             airTypeString,  1, 1, &file,    NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",             airTypeString,  1, 1, &ts,      NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output base name");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Convert HDF5 NIMROD Tokamak time step to NRRD format",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    tokamak_nimrod_parametric field(std::string(file), std::string(ts));
    
    std::ostringstream os;
    os << outs << "-t=" << ts << ".nrrd";
    
    Nrrd* nout = field.to_nrrd();
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell, nrrdCenterCell);
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                  << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    std::cout << "exported " << os.str() << std::endl;
}

}

