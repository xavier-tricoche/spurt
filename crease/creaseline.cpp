#include <vector>
#include <iostream>
#include <fstream>
#include <teem/hest.h>
#include <teem/limn.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <crease/measure_wrapper.hpp>
#include <crease/extractor.hpp>
#include <crease/pvo.hpp>
#include <crease/grid.hpp>

#include <math/math.hpp>
#include <image/probe.hpp>
#include <list>
#include <sstream>
#include <fstream>
#include <algorithm>

char* in_name, *flag_name, *out_name;
int measure, upsampling, maxdepth, maxdepthfix, ex_method;
bool ridge, trust_tracking = true;
float tolerance, grad_eps, v_thresh, v_thresh_select, s_thresh, s_thresh_select;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",                airTypeString,  1, 1, &in_name,             NULL,       "Input file name (NRRD format)");
    hestOptAdd(&hopt, "f",      "flags",                airTypeString,  1, 1, &flag_name,           NULL,       "Flag file name (NRRD format, 4D: value|strength)");
    hestOptAdd(&hopt, "o",      "output base",          airTypeString,  1, 1, &out_name,            NULL,       "Output file base name");
    hestOptAdd(&hopt, "m",      "measure",              airTypeInt,     0, 1, &measure,             "0",        "Invariant measure on data: 0: value, 1: FA, 2: mode");
    hestOptAdd(&hopt, "r",      "do ridge",             airTypeBool,    0, 1, &ridge,               "1",        "Extract ridges");
    hestOptAdd(&hopt, "d",      "max depth",            airTypeInt,     0, 1, &maxdepth,            "4",        "Max depth");
    hestOptAdd(&hopt, "D",      "max depth selecty",    airTypeInt,     0, 1, &maxdepthfix,         "6",        "Max depth in 2nd pass");
    hestOptAdd(&hopt, "u",      "upsampling",           airTypeInt,     0, 1, &upsampling,          "1",        "Spatial upsampling");
    hestOptAdd(&hopt, "e",      "extraction method",    airTypeInt,     0, 1, &ex_method,           "1",        "Extraction method", "0: PVO, 1: MC, 2: ZG");
    hestOptAdd(&hopt, "eps",    "tolerance",            airTypeFloat,   0, 1, &tolerance,           "0.05",     "Approximation tolerance");
    hestOptAdd(&hopt, "geps",   "gradient epsilon",     airTypeFloat,   0, 1, &grad_eps,            "1.0e-6",   "Min gradient norm for solution");
    hestOptAdd(&hopt, "v",      "value threshold",      airTypeFloat,   0, 1, &v_thresh,            "0",        "Value threshold for search");
    hestOptAdd(&hopt, "vs",     "solution threshold",   airTypeFloat,   0, 1, &v_thresh_select,     "0",        "Value threshold for solution");
    hestOptAdd(&hopt, "sv",     "strength threshold",   airTypeFloat,   0, 1, &s_thresh,            "0",        "Strength threshold for search");
    hestOptAdd(&hopt, "svs",    "strength solution",    airTypeFloat,   0, 1, &s_thresh_select,     "0",        "Strength threshold for solution");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize VCL curves in DPL dataset",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    using namespace spurt;
    
    initialize(argc, argv);
    
    spurt::crease::speedup = true;
    spurt::crease::read_info = true;
    spurt::crease::flag_file_name = flag_name;
    spurt::crease::bold_move = true;
    spurt::crease::display_debug_info = false;
    spurt::crease::crease_kind = measure;
    
    switch (ex_method) {
        case 0: {
            spurt::crease::ext_meth = spurt::crease::PVO;
            break;
        }
        case 1: {
            spurt::crease::ext_meth = spurt::crease::MC;
            break;
        }
        case 2: {
            spurt::crease::ext_meth = spurt::crease::ZG;
            break;
        }
    }
    
    spurt::crease::upsample = upsampling;
    
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, in_name, NULL)) {
        std::cerr << "unable to open " << in_name << std::endl;
        return -1;
    }
    crease::the_wrapper = new MeasureWrapper(nin, crease::crease_kind);
    
    std::vector< spurt::crease::line > lines;
    spurt::crease::vertices.clear();
    spurt::crease::value_threshold = v_thresh;
    spurt::crease::strength_threshold = s_thresh;
    spurt::crease::value_threshold_select = v_thresh_select;
    spurt::crease::strength_threshold_select = s_thresh_select;
    spurt::crease::confidence_threshold = 0.5;
    spurt::crease::max_depth = maxdepth;
    spurt::crease::max_depth_fix = maxdepthfix;
    spurt::crease::max_int_error = tolerance;
    spurt::crease::is_ridge = ridge;
    spurt::crease::upsample = upsampling;
    
    spurt::crease::extract_lines(lines, nin);
    
    // export results for reuse
    std::ostringstream os;
    os << out_name << "-" << (crease::is_ridge ? "ridge" : "valley")
       << "-up=" << upsampling
       << "-maxd=" << maxdepth << "-minval=" << crease::value_threshold
       << "-minstr=" << crease::strength_threshold << ".vcl";
    std::fstream output(os.str().c_str(), std::ios::out);
    if (output) {
        for (unsigned int l = 0 ; l < lines.size() ; l++) {
            for (unsigned int i = 0 ; i < lines[l].size() ; i++) {
                const vec3& p = spurt::crease::all_face_points[lines[l][i]];
                output << "p " << p[0] << " " << p[1] << " " << p[2] << std::endl;
            }
            output << "n" << std::endl;
        }
    }
    
    return 0;
}








