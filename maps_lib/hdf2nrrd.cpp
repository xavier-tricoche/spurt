#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <teem/nrrd.h>
#include <sstream>
#include <iostream>

char* fout, *fin, *ts;
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
    hestOptAdd(&hopt, "i",  "input",    airTypeString,  1, 1, &fin,     NULL,   "input hdf5 file name");
    hestOptAdd(&hopt, "t",  "time",     airTypeString,  1, 1, &ts,      NULL,   "time step");
    hestOptAdd(&hopt, "o",  "output",   airTypeString,  1, 1, &fout,    NULL,   "output nrrd file BASE name");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert the parametric representation of a HDF5 file to NRRD format",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(fin), std::string(ts));
    
    Nrrd* nrrd = field->to_nrrd();
    std::ostringstream os;
    os << fout << "-t=" << ts << ".nrrd";
    if (nrrdSave(os.str().c_str(), nrrd, NULL)) {
        std::cerr << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    return 0;
}
