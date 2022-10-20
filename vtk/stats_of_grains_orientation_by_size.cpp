#include <string>
#include <math/fixed_vector.hpp>
// #include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <teem/nrrd.h>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include <iomanip>
// #include "Garcia_vis_helper.hpp"

char*   param_in;
double  param_d;

// using namespace Garcia_vis_helper;


void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &param_in,            NULL,       "input size file (NRRD)");
    hestOptAdd(&hopt, "d",      "delta value",      airTypeDouble,  0, 1, &param_d,             "0.1",      "value increment (in %)");
    
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Visualize cyrstallographic orientation of granular microstructure",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    std::string size_name(param_in);
    
    Nrrd* __size = spurt::nrrd_utils::readNrrd(size_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> size(__size);
    std::cerr << size_name << " loaded.\n";
    int nb_grains = __size->axis[0].size;
    
    std::cerr << "there are " << nb_grains << " grains\n";
    
    // compute size range
    double min, max, d;
    std::vector<double> vals;
    {
        std::list<double> s;
        for (int i = 0 ; i < nb_grains ; ++i) {
            s.push_back(size[i]);
        }
        min = *std::min_element(s.begin(), s.end());
        max = *std::max_element(s.begin(), s.end());
        
        d = param_d * (max - min);
    }
    std::cerr << "min value = " << min << ", max value = " << max << ", delta val = " << d << '\n';
    
    int frame_counter = 0;
    for (double s = min ; s < max + 0.1*d ; s += d) {
        std::cerr << "\n\nFrame #" << frame_counter << ":\n";
        std::ostringstream os;
        os << "grains_orientation_on_sphere"
           << "_by_size_g=" << 1.0 << "_frame_" << std::setw(3) << std::setfill('0') << frame_counter++ << "_of_"
           << (int)floor((max - min) / d) + 1 << ".tiff";
           
        std::cerr << "exported frame: " << os.str() << '\n';
        
        double __min = s;
        double __max = s + d;
        
        std::cerr << "size range = [" << __min << ", " << __max << "]\n";
        
        // included grains
        int local_min, local_max;
        float local_avg;
        local_min = std::numeric_limits<int>::max();
        local_max = 0;
        local_avg = 0;
        std::set<int> selected_grains;
        for (int i = 0 ; i < nb_grains ; ++i) {
            int sz = size[i];
            if (sz >= __min && sz <= __max) {
                selected_grains.insert(i);
                if (sz < local_min) {
                    local_min = sz;
                } else if (sz > local_max) {
                    local_max = sz;
                }
                local_avg += sz;
            }
        }
        
        int N = selected_grains.size();
        
        std::cerr << "statistics:\n"
                  << "\t number of selected grains: " << N << " (" << N*100 / nb_grains << "%)\n"
                  << "\t min grain size: " << local_min << '\n'
                  << "\t max grain size: " << local_max << '\n'
                  << "\t average grain size: " << local_avg / (float)N << '\n';
    }
}















































