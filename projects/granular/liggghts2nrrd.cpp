#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <teem/nrrd.h>
#include <vector>

char* in_name, *out_base;
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
    hestOptAdd(&hopt, "i",      "input file",       airTypeString,  1, 1, &in_name,     NULL,   "input file name (LIGGGHTS dump file)");
    hestOptAdd(&hopt, "o",      "output base",      airTypeString,  1, 1, &out_base,    NULL,   "output base name");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert LIGGGHTS dump file to NRRD format",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

bool read_timestep(std::fstream& f, size_t& ti, size_t& np, std::vector<std::string>& attributes, float** data)
{
    static const std::string ts_str("ITEM: TIMESTEP");
    static const std::string na_str("ITEM: NUMBER OF ATOMS");
    static const std::string bb_str("ITEM: BOX BOUNDS");
    static const std::string at_str("ITEM: ATOMS");
    std::string buffer, dummy, var;
    std::getline(f, buffer); // ITEM: TIMESTEP\n
    if (buffer.compare(0, ts_str.size(), ts_str)) {
        return false;
    }
    f >> ti;
    std::getline(f, buffer); // '\n'
    std::getline(f, buffer); // ITEM: NUMBER OF ATOMS\n
    if (buffer.compare(0, na_str.size(), na_str)) {
        return false;
    }
    f >> np;
    std::getline(f, buffer); // '\n'
    std::getline(f, buffer); // ITEM: BOX BOUNDS pp pp ff\n
    if (buffer.compare(0, bb_str.size(), bb_str)) {
        return false;
    }
    std::getline(f, buffer); // xlo xhi\n
    std::getline(f, buffer); // ylo yhi\n
    std::getline(f, buffer); // zlo zhi\n
    std::getline(f, buffer); // ITEM: ATOMS ...\n
    if (buffer.compare(0, at_str.size(), at_str)) {
        return false;
    }
    std::istringstream is(buffer);
    is >> dummy >> dummy >> var; // ITEM: ATOMS
    attributes.clear();
    while (!is.eof() && var !="\n") {
        attributes.push_back(var);
        is >> var;
    }
    (*data) = (float*)calloc(attributes.size()*np, sizeof(float));
    for (int i=0 ; i<attributes.size()*np ; ++i) {
        f >> (*data)[i];
    }
    std::getline(f, buffer); // '\n'
    return true;
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    std::fstream input(in_name, std::ios::in);
    if (!input) {
        std::cerr << argv[0] << ": unable to open " << in_name;
        exit(-1);
    }
    
    std::ostringstream os;
    float* data;
    for (int c=0; !input.eof() ; ++c) {
        os.clear();
        os.str("");
        os << out_base << "_" << std::setw(6) << std::setfill('0') << c << ".nrrd";
        size_t ti, np;
        std::vector<std::string> attributes;
        if (read_timestep(input, ti, np, attributes, &data)) {
            Nrrd* nout = nrrdNew();
            size_t size[] = {attributes.size(), np};
            if (nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, size)) {
                std::cerr << biffGetDone(NRRD) << std::endl;
                exit(-1);
            }
            std::ostringstream os2;
            const char* labels[2];
            for (int i=0 ; i<attributes.size() ; ++i) {
                os2 << attributes[i].c_str();
                if (i<attributes.size()-1) {
                    os2 << ';';
                }
            }
            labels[0] = os2.str().c_str();
            labels[1] = "particles";
            nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
            if (nrrdSave(os.str().c_str(), nout, NULL)) {
                std::cerr << biffGetDone(NRRD) << std::endl;
                exit(-1);
            }
            std::cerr << "exported:\t" << os.str() << std::endl;
            nrrdNuke(nout);
        } else {
            std::cerr << "time step prematurily interrupted. exit" << std::endl;
            exit(-1);
        }
    }
    
    return 0;
    
}


















