#include <vector>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <image/nrrd_wrapper.hpp>
#include <misc/progress.hpp>
#include <math/basic_math.hpp>
#include <math/small_tensor.hpp>
#include <math/types.hpp>

char* outs, *file, *att;
bool verbose;

using namespace spurt;

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
    hestOptAdd(&hopt, "i",      "file",             airTypeString,  1, 1, &file,    NULL,       "input file name (NRRD)");
    hestOptAdd(&hopt, "a",      "attribute file",   airTypeString,  0, 1, &att,     "none",     "attribute file name (NRRD)");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output file base name");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &verbose, "0",        "verbose mode (debugging)");
    
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute statistical properties of a granular microstructure",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline ivec3 int_to_ivec(int i, const svec3& s)
{
    int x = i % s[0];
    int aux = i / s[0];
    int y = aux % s[1];
    int z = aux / s[1];
    return ivec3(x, y, z);
}

inline vec3 coord(int i, const svec3& size, const vec3& step)
{
    vec3 c = int_to_ivec(i, size);
    return c * step;
}

int main(int argc, const char* argv[])
{
	initialize(argc, argv);
	std::string basename(outs);

	Nrrd* nin = spurt::nrrd_utils::readNrrd(file);
	std::vector<int> tags;
	spurt::nrrd_utils::to_vector(tags, nin);
	svec3 size;
	vec3 step;
	for (int i = 0 ; i < 3 ; ++i) {
		size[i] = nin->axis[i].size;
		step[i] = nin->axis[i].spacing;
		if (!step[i] || std::isnan(step[i]) || std::isinf(step[i])) step[i] = 1;
	}
	std::cout << "step = " << step << '\n';
	std::cout << "size = " << size << '\n';

	bool use_attributes = false;
	std::vector<float> attributes;
	if (strcmp(att, "none")) {
		use_attributes = true;
		Nrrd *tmp = spurt::nrrd_utils::readNrrd(att);
		spurt::nrrd_utils::to_vector(attributes, tmp);
	}

	// basic objects of the processing
	// 1. voxels of the input dataset
	// 2. grain tags (value associated with each voxel)
	// 3. grain ids (id associated with each grain in output)

	size_t nbvoxels = size[0] * size[1] * size[2];

	typedef	int										Key;
	typedef std::pair<Key, int>						value_type;
	typedef std::multimap<Key, int>					multimap_type;
	typedef multimap_type::iterator					iterator_type;
	typedef std::pair<iterator_type, iterator_type>	range_type;

	// invert the relationship voxel to tag
	multimap_type tag_to_voxel;
	// assign ids to grain tags
	std::map<int, int> tag_to_id;

	size_t nb_tags = 0;
	for (int i = 0 ; i < tags.size() ; ++i) {
		int tag = tags[i];
		if (tag_to_id.find(tag) == tag_to_id.end()) {
			tag_to_id[tag] = nb_tags++;
		}
	}
	std::cerr << "there are " << nb_tags << " grains in input\n";
	{
		std::string name(basename);
		name.append("-ids_to_tags.nrrd");
		Nrrd *nout = nrrdNew();
		int *ids = (int*)calloc(nb_tags, sizeof(int));
		for (std::map<int, int>::iterator i = tag_to_id.begin() ; i != tag_to_id.end() ; ++i) {
			ids[i->second] = i->first;
		}

		size_t sz = nb_tags;
		if (nrrdWrap_nva(nout, ids, nrrdTypeInt, 1, &sz)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	for (int i = 0 ; i < tags.size() ; ++i) {
		tag_to_voxel.insert(value_type(tags[i], i));
	}

	// compute grain sizes
	int *grain_size = (int*)calloc(nb_tags, sizeof(int));
	for (std::map<int, int>::iterator i = tag_to_id.begin() ; i != tag_to_id.end() ; ++i) {
		grain_size[i->second] = tag_to_voxel.count(i->first);
	}


	int *grain_id_to_size = (int*)calloc(2 * nb_tags, sizeof(int));
	{
		int count = 0;
		for (multimap_type::iterator i = tag_to_voxel.begin() ;
		i != tag_to_voxel.end() ; i = tag_to_voxel.upper_bound(i->first), ++count) {
			grain_id_to_size[2*count] = i->first;
			grain_id_to_size[2*count+1] = tag_to_voxel.count(i->first);
		}
	}

	if (use_attributes) {
		float *minmaxs = (float*)calloc(2 * nb_tags, sizeof(float));
		for (int i = 0 ; i < nb_tags ; ++i) {
			minmaxs[2*i] = std::numeric_limits<float>::max();
			minmaxs[2*i+1] = std::numeric_limits<float>::min();
		}

		for (int i = 0 ; i < nbvoxels ; ++i) {
			int id = tag_to_id[tags[i]];
			float val = attributes[i];
			minmaxs[2*id] = std::min(minmaxs[2*id], val);
			minmaxs[2*id+1] = std::max(minmaxs[2*id+1], val);
		}

		std::string name(basename);
		name.append("-grain-span.nrrd");
		Nrrd *nout = nrrdNew();
		size_t _size[] = {2, nb_tags};
		if (nrrdWrap_nva(nout, minmaxs, nrrdTypeFloat, 2, _size)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	// debug
	{
		std::string name(basename);
		name.append("-size-inplace.nrrd");
		Nrrd *nout = nrrdNew();
		size_t *__sizes = (size_t*)calloc(tags.size(), sizeof(size_t));
		for (int i = 0 ; i < tags.size() ; ++i) {
			__sizes[i] = grain_size[tag_to_id[tags[i]]];
		}

		size_t _size[3] = { size[0], size[1], size[2] };
		if (nrrdWrap_nva(nout, __sizes, nrrdTypeInt, 3, _size)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	{
		std::string name(basename);
		name.append("-size2.nrrd");
		Nrrd *nout = nrrdNew();
		size_t sz = nb_tags;
		if (nrrdWrap_nva(nout, grain_size, nrrdTypeInt, 1, &sz)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	{
		std::string name(basename);
		name.append("-size.nrrd");
		Nrrd *nout = nrrdNew();
		size_t sz[2] = {2, nb_tags};
		if (nrrdWrap_nva(nout, grain_id_to_size, nrrdTypeInt, 2, sz)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	// compute grain anisotropy
	float *voxel_fa = (float*)calloc(nb_tags, sizeof(float));
	float *voxel_aniso = (float*)calloc(nb_tags * 3, sizeof(float));
	float *voxel_dir = (float*)calloc(nb_tags * 3, sizeof(float));
	for (iterator_type it = tag_to_voxel.begin() ; it != tag_to_voxel.end() ;
	     it = tag_to_voxel.upper_bound(it->first)) {
		int tag = it->first;
		range_type range = tag_to_voxel.equal_range(tag);
		vec3 center(0, 0, 0);
		for (iterator_type i = range.first ; i != range.second ; ++i) {
			vec3 x = coord(i->second, size, step);
			center += x;
		}
		center /= (float)tag_to_voxel.count(tag);

		spurt::mat3 A;
		for (iterator_type i = range.first ; i != range.second ; ++i) {
			vec3 x = coord(i->second, size, step) - center;
			A += spurt::outer(x,x);
		}
		int id = tag_to_id[tag];
		voxel_fa[id] = spurt::FA(A);
		if (std::isnan(voxel_fa[id]) || std::isinf(voxel_fa[id])) voxel_fa[id] = 0;
		vec3 aniso = spurt::westin_aniso(A);
		for (int j = 0 ; j < 3 ; ++j) voxel_aniso[3*id+j] = aniso[j];

        cvec3 evals;
        cmat3 evecs;
        spurt::eigensystem(evals, evecs, A);
        vec3 emax = spurt::real(evecs.column(0));
		for (int j = 0 ; j < 3 ; ++j) voxel_dir[3*id+j] = emax[j];
	}

	// debug
	{
		std::string name(basename);
		name.append("-fa-inplace.nrrd");
		Nrrd *nout = nrrdNew();
		float *__fa = (float*)calloc(tags.size(), sizeof(float));
		for (int i = 0 ; i < tags.size() ; ++i) {
			__fa[i] = voxel_fa[tag_to_id[tags[i]]];
		}

		size_t _size[3] = { size[0], size[1], size[2] };
		if (nrrdWrap_nva(nout, __fa, nrrdTypeFloat, 3, _size)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	{
		std::string name(basename);
		name.append("-fa.nrrd");
		Nrrd *nout = nrrdNew();
		size_t sz = nb_tags;
		if (nrrdWrap_nva(nout, voxel_fa, nrrdTypeFloat, 1, &sz)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	// debug
	{
		std::string name(basename);
		name.append("-westin-inplace.nrrd");
		Nrrd *nout = nrrdNew();
		float *__aniso = (float*)calloc(3 * tags.size(), sizeof(float));
		for (int i = 0 ; i < tags.size() ; ++i) {
			int id = tag_to_id[tags[i]];
			__aniso[3*i  ] = voxel_aniso[3*id  ];
			__aniso[3*i+1] = voxel_aniso[3*id+1];
			__aniso[3*i+2] = voxel_aniso[3*id+2];
		}

		size_t _size[4] = { 3, size[0], size[1], size[2] };
		if (nrrdWrap_nva(nout, __aniso, nrrdTypeFloat, 4, _size)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}
	{
		std::string name(basename);
		name.append("-westin.nrrd");
		Nrrd *nout = nrrdNew();
		size_t sz[] = {3, nb_tags};
		if (nrrdWrap_nva(nout, voxel_aniso, nrrdTypeFloat, 2, sz)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}
	{
		std::string name(basename);
		name.append("-direction.nrrd");
		Nrrd *nout = nrrdNew();
		size_t sz[] = {3, nb_tags};
		if (nrrdWrap_nva(nout, voxel_dir, nrrdTypeFloat, 2, sz)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		if (nrrdSave(name.c_str(), nout, NULL)) {
			std::cerr << "grain-stat: " << biffGetDone(NRRD) << std::endl;
			exit(-1);
		}
		nrrdNuke(nout);
	}

	return 0;
}























