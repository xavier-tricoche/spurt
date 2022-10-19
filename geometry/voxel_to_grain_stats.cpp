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

char* outs, *file, *att, *atts;
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
    hestOptAdd(&hopt, "i",      "file",             airTypeString,  1, 1, &file,    NULL,       "input file name (NRRD)");
    hestOptAdd(&hopt, "a",      "attribute file",   airTypeString,  1, 1, &att,     NULL,       "attribute file name (NRRD)");
    hestOptAdd(&hopt, "as",     "attribute name",   airTypeString,  1, 1, &atts,    NULL,       "attribute description");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output file base name");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &verbose, "0",        "verbose mode (debugging)");
    
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute statistical properties of a granular microstructure",
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
    std::string basename(outs);
    
    Nrrd* nin = xavier::readNrrd(file);
    std::vector<int> tags;
    xavier::to_vector(tags, nin);
    nvis::ivec3 size;
    nvis::vec3 step;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i] = nin->axis[i].size;
        step[i] = nin->axis[i].spacing;
        if (!step[i] || std::isnan(step[i]) || std::isinf(step[i])) {
            step[i] = 1;
        }
    }
    std::cout << "step = " << step << '\n';
    std::cout << "size = " << size << '\n';
    
    std::vector<float> attributes;
    Nrrd* tmp = xavier::readNrrd(att);
    xavier::to_vector(attributes, tmp);
    
    int nbvoxels = size[0] * size[1] * size[2];
    
    typedef int                                     Key;
    typedef std::pair<Key, int>                     value_type;
    typedef std::multimap<Key, int>                 multimap_type;
    typedef multimap_type::iterator                 iterator_type;
    typedef std::pair<iterator_type, iterator_type> range_type;
    
    // invert the relationship voxel to tag
    multimap_type tag_to_voxel;
    
    for (int i = 0 ; i < nbvoxels ; ++i) {
        tag_to_voxel.insert(value_type(tags[i], i));
    }
    // compute grain sizes
    int nb_tags = 0;
    for (iterator_type it = tag_to_voxel.begin() ; it != tag_to_voxel.end() ;
            it = tag_to_voxel.upper_bound(it->first)) {
        ++nb_tags;
    }
    int* grain_size = (int*)calloc(2 * nb_tags, sizeof(int));
    int counter = 0;
    for (iterator_type it = tag_to_voxel.begin() ; it != tag_to_voxel.end() ;
            it = tag_to_voxel.upper_bound(it->first)) {
        int id = it->first;
        int sz = tag_to_voxel.count(id);
        grain_size[counter++] = id;
        grain_size[counter++] = sz;
    }
    Nrrd* nout = nrrdNew();
    size_t s1[] = {2, nb_tags};
    nrrdWrap_nva(nout, grain_size, nrrdTypeInt, 2, s1);
    std::string name(basename);
    name.append("-grain-size.nrrd");
    nrrdSave(name.c_str(), nout, NULL);
    nrrdNuke(nout);
    
    // value_type = <ID, min, max>
    float* minmaxs = (float*)calloc(3 * nb_tags, sizeof(float));
    float* averages = (float*)calloc(2 * nb_tags, sizeof(float));
    std::map<int, std::pair<float, float> > grain_span;
    std::map<int, std::pair<float, int> > grain_sum;
    for (iterator_type it = tag_to_voxel.begin() ; it != tag_to_voxel.end() ; ++it) {
        int id = it->first;
        float val = attributes[it->second];
        if (grain_span.find(id) == grain_span.end()) {
            grain_span[id] = std::pair<float, float>(val, val);
            grain_sum[id] = std::pair<float, int>(val, 1);
        } else {
            std::pair<float, float>& _span = grain_span[id];
            if (val < _span.first) {
                _span.first = val;
            } else if (val > _span.second) {
                _span.second = val;
            }
            std::pair<float, int>& _sum = grain_sum[id];
            _sum.first += val;
            ++_sum.second;
        }
    }
    counter = 0;
    for (std::map<int, std::pair<float, float> >::const_iterator it = grain_span.begin();
            it != grain_span.end() ; ++it) {
        minmaxs[counter++] = it->first;
        minmaxs[counter++] = it->second.first;
        minmaxs[counter++] = it->second.second;
    }
    
    size_t s2[] = {3, nb_tags};
    nout = nrrdNew();
    nrrdWrap_nva(nout, minmaxs, nrrdTypeFloat, 2, s2);
    name = basename;
    std::ostringstream os;
    os << basename << "-grain-" << atts << "-span.nrrd";
    nrrdSave(os.str().c_str(), nout, NULL);
    nrrdNuke(nout);
    
    counter = 0;
    for (std::map<int, std::pair<float, int> >::const_iterator it = grain_sum.begin();
            it != grain_sum.end() ; ++it) {
        averages[counter++] = it->first;
        averages[counter++] = it->second.first / (float)it->second.second;
    }
    
    s2[0] = 2;
    nout = nrrdNew();
    nrrdWrap_nva(nout, averages, nrrdTypeFloat, 2, s2);
    name = basename;
    os.clear();
    os.str("");
    os << basename << "-grain-" << atts << "-average.nrrd";
    nrrdSave(os.str().c_str(), nout, NULL);
    nrrdNuke(nout);
    
    return 0;
}


































