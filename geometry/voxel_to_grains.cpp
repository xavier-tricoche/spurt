#include <vector>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>
#include <util/timer.hpp>
#include <math/math.hpp>
#include <math/matrix.hpp>
#include <math/tensor.hpp>

char* outs, *tags, *vals;
bool verbose;
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
    hestOptAdd(&hopt, "t",      "tags",         airTypeString,  1, 1, &tags,    NULL,       "input tag file (NRRD)");
    hestOptAdd(&hopt, "v",      "values",       airTypeString,  1, 1, &vals,    NULL,       "input value file (NRRD)");
    hestOptAdd(&hopt, "o",      "output",       airTypeString,  1, 1, &outs,    NULL,       "output file name");
    
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Turn voxel-based into grain-based values",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline nvis::ivec3 int_to_ivec(int i, const nvis::ivec3& s)
{
    int x = i % s[0];
    int aux = i / s[0];
    int y = aux % s[1];
    int z = aux / s[1];
    return nvis::ivec3(x, y, z);
}

inline nvis::vec3 coord(int i, const nvis::ivec3& size, const nvis::vec3& step)
{
    nvis::vec3 c = int_to_ivec(i, size);
    return c*step;
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin = spurt::readNrrd(tags);
    std::vector<int> __tags;
    spurt::to_vector(__tags, nin);
    
    std::vector<float> __vals;
    Nrrd* tmp = spurt::readNrrd(vals);
    spurt::to_vector(__vals, tmp);
    
    int nbvoxels = __tags.size();
    
    std::map<int, nvis::vec3>   grain_values;
    
    for (int i = 0 ; i < nbvoxels ; ++i) {
        int id = __tags[i];
        if (grain_values.find(id) == grain_values.end()) {
            grain_values[id] = nvis::vec3(__vals[3*i], __vals[3*i+1], __vals[3*i+2]);
        }
    }
    
    int nb_grains = grain_values.size();
    float* out = (float*)calloc(4 * nb_grains, sizeof(float));
    int id = 0;
    for (std::map<int, nvis::vec3>::const_iterator it = grain_values.begin() ;
            it != grain_values.end() ; ++it) {
        int k = it->first;
        const nvis::vec3& v = it->second;
        out[id++] = k;
        for (int i = 0 ; i < 3 ; ++i) {
            out[id++] = v[i];
        }
    }
    
    size_t sz[] = {4, nb_grains};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, out, nrrdTypeFloat, 2, sz);
    nrrdSave(outs, nout, NULL);
    
    return 0;
}

























