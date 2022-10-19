#include <iostream>
#include <sstream>

#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/wall_timer.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;

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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1,  1,  &name_in,       NULL,   "input file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Export a RGB color map of a FTLE NRRD",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef xavier::nrrd_data_traits<Nrrd*>  field_type;

inline nvis::vec3 color(double f, double b, double min, double max)
{
    static nvis::vec3 red(1,0,0);
    static nvis::vec3 white(1,1,1);
    static nvis::vec3 blue(0,0,1);
    static nvis::vec3 black(0,0,0);
    
    double u = (f-min)/(max-min);
    double v = (b-min)/(max-min);
    
    return (1.0-u)*(1.0-v)*white + u*(1.0-v)*red + u*v*black + (1.0-u)*v*blue;
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin = xavier::nrrd_utils::readNrrd(name_in);
    std::vector<size_t> size(nin->dim);
    for (int i=0 ; i<size.size() ; ++i) {
        size[i] = nin->axis[i].size;
    }
    int nb_channels = size[0];
    int nbpix = 1;
    for (int i=1 ; i<size.size() ; ++i) {
        nbpix *= size[i];
    }
    
    double max;
    std::vector<float> vals;
    xavier::nrrd_utils::to_vector<float>(vals, nin);
    max = *std::max_element(vals.begin(), vals.end());
    std::cerr << "max = " << max << std::endl;
    
    int shift = nb_channels-1;
    
    float* out = (float*)calloc(nbpix*3, sizeof(float));
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < nbpix ; ++n) {
            nvis::fvec3 c = color(vals[nb_channels*n], vals[nb_channels*n+shift], 0, max);
            out[3*n  ] = c[0];
            out[3*n+1] = c[1];
            out[3*n+2] = c[2];
        }
    }
    
    std::cout << "done" << std::endl;
    
    size[0] = 3;
    std::vector<double> spc(size.size());
    spc[0] = airNaN();
    for (int i=1 ; i<size.size() ; ++i) {
        spc[i] = nin->axis[i].spacing;
    }
    xavier::nrrd_utils::writeNrrdFromContainers(out, name_out, /*nrrdTypeFloat,*/ size, spc);
    
    return 0;
}


























