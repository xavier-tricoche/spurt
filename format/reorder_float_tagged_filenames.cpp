#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <teem/nrrd.h>
#include <map>
#include <iomanip>

char* input;
char* prefix;
char* suffix;
char* output;

void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",    airTypeString,  1,  1,  &input,     NULL,   "input file name");
    hestOptAdd(&hopt, "o",      "output",   airTypeString,  1,  1,  &output,    NULL,   "output script");
    hestOptAdd(&hopt, "pre",    "prefix",   airTypeString,  1,  1,  &prefix,    NULL,   "prefix");
    hestOptAdd(&hopt, "suf",    "suffix",   airTypeString,  1,  1,  &suffix,    NULL,   "sufffix");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert file names",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, const char* argv[])
{
    initialize(argc, argv);
    
    std::string pre(prefix), suf(suffix);
    
    std::map<float, std::string> tag_to_name;
    
    std::fstream in(input, std::ios::in);
    std::fstream out(output, std::ios::out);
    
    out << "#!/bin/sh\n";
    while (!in.eof()) {
        std::string name;
        in >> name;
        
        int starts_at = pre.size();
        int ends_at = name.rfind(suf);
        std::string val_str;
        try {
            val_str.assign(name, starts_at, ends_at - starts_at);
        } catch (...) {
            break;
        }
        float tag = atof(val_str.c_str());
        std::cerr << name << " -> " << val_str << " -> " << tag << '\n';
        tag_to_name[tag] = name;
    }
    
    std::cerr << "re-ordered filenames are:\n";
    int id = 0;
    for (std::map<float, std::string>::iterator it = tag_to_name.begin() ; it != tag_to_name.end() ; ++it) {
        std::cerr << it->second << '\n';
        
        out << "mv " << it->second << " img=" << id++ << "_" << it->second << '\n';
    }
    out.close();
    in.close();
    
    return 0;
}

