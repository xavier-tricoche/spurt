#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include <math/fixed_vector.hpp>

#include <granular/dem.hpp>
#include <granular/dem_utils.hpp>
#include <granular/boundary.hpp>
#include <granular/neighbors.hpp>
#include <granular/particles.hpp>
#include <granular/walton_forces.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr nvis::vec3 gravity(0, 0, -9.80665);

using namespace xavier::granular::dem;

/******************************************************************************
*                                                                             *
*                              TYPE DEFINITIONS                               *
*                                                                             *
******************************************************************************/                        
// basic types
typedef size_t                                     index_t;
typedef double                                     value_t;
typedef nvis::fixed_vector<value_t, 3>             vector_t;
typedef nvis::fixed_vector<bool, 3>                bvec3;
typedef nvis::bounding_box<vector_t>               bounds_t;
typedef std::pair<index, index>                    id_pair;
typedef std::pair<vector_t, vector_t>              force_pair;

// particle-related
typedef Spherical<index, value_t, vector_t>        particle;
typedef ParticleTraits<particle>                   particle_traits;

// boundary-related
typedef PlanarBoundary<index_t, value_t, vector_t> boundary;
typedef BoundaryTraits<boundary>                   boundary_traits;

// force-related
typedef walton::ContactInformation<index_t, value_t, vector_t> 
                                                   contact_information;
typedef walton::ContactHistory<vector_t, value_t>  history;
typedef walton::ForceModel<index_t, value_t, vector_t, particle, 
                           particle_traits, boundary, boundary_traits,
                           history>                force_model;

// data structure
typedef CellGrid<index, particle>                  cell_grid;
typedef NeighborList<index, particle>              neighbor_list;
typedef std::vector<particle>                      particle_container;
typedef std::vector<neighbor_list>                 neighbor_container;
typedef std::map<contact, force_pair>              force_container;
typedef std::map<contact, history>                 history_container;
typedef force_container::const_iterator            const_force_iterator;
typedef force_container::iterator                  force_iterator;
typedef history_container::const_iterator          const_history_iterator;
typedef history_container::iterator                history_iterator;

// particle information
particle_container __particles;
neighbor_container __neighbor_lists;

// tap parameters
value_t  __a     = 0.01; // tap amplitude
value_t  __f     = 7.5;  // tap frequency
value_t  __ntaps = 5;    // number of taps
value_t  __T     = 1;    // tapping period

// domain parameters
bounds_t __bounds;                              // domain boundary
vector_t __bmin     = vector_t(0, 0, 0);        // low bound
vector_t __bmax     = vector_t(0.24, 0.24, 10); // high bound
bvec3    __periodic = bvec3(true, true, false); // domain periodicity

// force model
value_t  __gamma    = 1./3.;  // Mindlin's elastic friction parameter
value_t  __k_load   = 280000; // loading spring constant
value_t  __k_unload = 330810; // unloading spring constant
value_t  __k_ratio  = 0.8;    // tangential vs. normal stiffness
value_t  __mu       = 0.1;    // coefficient of friction

// particle properties
value_t  __density = 1200; // particle density
value_t  __mass;           // particle mass
value_t  __radius  = 0.01; // particle radius

// spatial data structure parameters
value_t  __width = 0.03;  // cell width in cell grid
value_t  __skin  = 0.005; // skin thickness

// computation parameters
int      __nthreads = 0;      // maximum available
value_t  __dt       = 1.0e-7; // integration step size

// I/O
std::string          __input_name;
std::string          __output_basename;
std::string          __log_name         = "-"; // log output to stdio
null_ostream_buffer* __null_buffer_ptr  = 0;
dual_ostream_buffer* __dual_buffer_ptr  = 0;
dual_sink*           __dual_sink_ptr    = 0;
std::streambuf*      __cerr_buffer_ptr  = 0;
std::ostream         __log;
std::string          __me;
float                __output_freq      = 100;
bool                 __verbose          = false;

/*****************************************************************************
*                                                                            *
*                                I/O FUNCTIONS                               *
*                                                                            *
*****************************************************************************/                        

void clean_exit(int exit_code=0) {
    // garbage collection
    if (__null_buffer_ptr) delete __null_buffer_ptr;
    if (__dual_buffer_ptr) delete __dual_buffer_ptr;
    if (__dual_sink_ptr)   delete __dual_sink_ptr;
    // restore std::cerr's stream buffer
    if (__cerr_buffer_ptr) std::cerr.rdbuf(__cerr_buffer_ptr);
    
    exit(exit_code);
}

void warning(const std::string& what, bool do_exit=true) {
    std::cerr << "\n\nERROR: " << what << std::endl;
    if (do_exit) clean_exit(-1);
}

void usage(const std::string& msg="")
{
    if (!msg.empty()) warning("Missing " + what, false);
    
    __log
        << "USAGE: " << me << " <parameters> [<options>]\n"
        << "PARAMETERS:\n"
        << " -i  | --input <string>      Input file name\n"
        << " -o  | --output <string>     Output base name\n"
        << "OPTIONS:\n"                  
        << " -h  | --help                Print this information\n"
        << " -v  | --verbose             Verbose mode\n"
        << " \% Particle Properties...\n"
        << " -r  | --radius <float>      Particle radius (default: 0.01 m)\n"
        << " -d  | --density <float>     Material density (default: 1200 kg/m^3)\n"
        << " \% Domain Properties...\n"
        << " -b  | --boundary <float>x6  Simulation domain boundaries (default: 0 0.24 0 0.24 0 1)\n"
        << " -p  | --periodic <bool>x3   Domain periodicity (default: true true false)\n"
        << " \% Force Model Parameters...\n"
        << " -kn | --knormal <float>x2   Loading/unloading spring stiffness (default: 28000/33081 N/m)\n"
        << " -kr | --kratio <float>      Tangential vs. normal stiffness ratio (default: 0.8)\n"
        << " -g  | --gamma <float>       Mindlin's elastic friction parameter (default: 0.333...)\n"
        << " -m  | --mu <float>          Coefficient of friction (default: 0.1)\n"
        << " \% Tapping Parameters...\n"
        << " -f  | --frequency <float>   Tap frequency (default: 7.5 Hz)\n"
        << " -a  | --amplitude <float>   Tap amplitude (default: 0.02 m)\n"
        << " -T  | --period <float>      Time interval between taps (default: 1 s)\n"
        << " -n  | --number <int>        Number of taps (default: 5)\n"
        << " \% Computation Parameters...\n"
        << " -w  | --width <float>       Cell width in cell grid (default: 0.025 m)\n"
        << " -s  | --skin <float>        Skin thickness (default: 0.005 m)\n"
        << " -dt | --step <float>        Integration step size (default: 1e-7 s)\n"
        << " \% I/O Parameter...\n"
        << " -l  | --log <string>        Log / diagnostics output file name (default: stdout \"-\")\n"
        << " -of | --outputf <float>     Output frequency (default: 100 Hz)\n"
        << std::endl;
    
    clean_exit(what.empty() ? 0 : -1);
}

namespace io = namespace boost::iostreams;
typedef io::basic_null_sink<char>                   null_sink;
typedef io::stream_buffer<nullsink>                 null_ostream_buffer;
typedef io::tee_device<std::ostream, std::ofstream> dual_sink;
typedef io::stream_buffer<dual_sink>                dual_ostream_buffer;

template<typename Target, typename Source>
struct lexical_cast_nothrow {
    static bool operator(const Source& s) const {
        try {
            Target cast_ = boost::lexical_cast<Target>(s);
        }
        catch(boost::bad_lexical_cast& e) {
            return false;
        }
        return true;
    }
};
typedef lexical_cast_nothrow<int, std::string>    is_an_int;
typedef lexical_cast_nothrow<double, std::string> is_a_float;

bool import_data_plain(std::fstream& input) {
    index_t nb_particles = 0;
    std::string end_of_line, keyword;
    char c;
    
    // read initial configuration lines
    input.peek(c);
    while (c == '#') {
        input >> c; // '#'
        input >> keyword;
        if (keyword == "number") {
            input >> nb_particles;
            if (__verbose) 
        }
        else if (keyword == "X-range") {
            input >> __bounds.min()[0] >> __bounds.max()[0];
        }
        else if (keyword == "Y-range") {
            input >> __bounds.min()[1] >> __bounds.max()[1];
        }
        else if (keyword == "Z-range") {
            input >> __bounds.min()[2] >> __bounds.max()[2];
        }
        else if (keyword == "radius") {
            input >> __radius;
        }
        else if (keyword == "periodic") {
            input >> __periodic[0] >> __periodic[1] >> __periodic[2];
        }
        input.getline(end_of_line); // skip '\n'
        input.peek(c);
    }
    
    __particles.resize(nb_particles);
    for (index_t i=0 ; i<nb_particles ; ++i) {
        particle& p = __particles[i];
        p.index = i;
        input >> p.position[0] >> p.position[1] >> p.position[2];
        input.getline(end_of_line);
    }
    input.close();
    
    return true;
}

bool import_data_liggghts(std::fstream& input) {
    index_t nb_particles=0, nb_atom_types=1;
    std::string end_of_line, keyword;
    
    // read initial configuration lines
    while (!input.eof()) {
        std::string line, word;
        std::vector<std::string> words;
        
        input.getline(line);        
        // skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        const std::string error_msg = 
            "Trouble reading input file.\n" \
            "Invalid format in line:\n\t\"" + line + "\"";
        
        // otherwise parse line
        std::istringstream iss(line);
        while (!iss.eof()) {
            iss >> word;
            words.push_back(word);
        }
        words.pop_back();
        
        // check line contents
        // "<N> atoms"
        if (words.size() == 2 && words[1] == "atoms") {
            nb_particles = atoi(words[0].c_str());
        }
        // "<T> atom types"
        else if (words.size() == 3 &&
                 is_an_int(words[0]) && 
                 words[1] == "atom" && words[2] == "types") {
            nb_atom_types = atoi(words[0]);
        }
        else if (words.size() == 4 && 
                 is_a_float(words[0]) && is_a_float(words[1])) {
            // "<xmin> <xmax> xlo xhi"
            if (words[2] == "xlo" && words[3] == "xhi") {
                __bounds.min()[0] = atof(words[0]);
                __bounds.max()[0] = atof(words[1]);
            }
            // "<ymin> <ymax> ylo yhi"
            else if (words[2] == "ylo" && words[1] == "yhi") {
                __bounds.min()[1] = atof(words[0]);
                __bounds.max()[1] = atof(words[1]);
            }
            // "<zmin> <zmax> zlo zhi"
            else if (words[2] == "zlo" && words[3] == "zhi") {
                __bounds.min()[2] = atof(words[0]);
                __bounds.max()[2] = atof(words[1]);
            }
            else warning(error_msg);
        }
        // "Atoms"
        else if (words.size() == 1 && words[0] == "Atoms") {
            // skip arbitrary number of empty lines
            input.getline(line);
            while(line.empty()) { input.getline(); }
            break;
        }
        else warning(error_msg);
    }

    // read particle information, assumed constant for all particles
    __particles.resize(nb_particles);
    index_t _type;
    value_t _rad, _dens;
    input >> __particles[0].index >> _type >> _rad >> _dens 
          >> __particles[0].position[0]
          >> __particles[0].position[1]
          >> __particles[0].position[2];
    input.getline(end_of_line);
    
    __radius = _rad;
    __density = _dens;
    
    for (index_t i=0 ; i<nb_particles ; ++i) {
        std::string line, word;
        std::vector<std::string> words;
        particle& p = __particles[i];
        input >> p.index >> _type >> _rad >> dens
              >> p.position[0] >> p.position[1] >> p.position[2];
        input.getline(end_of_line);
    }
    input.close();
    
    return true;
}

int main(int argc, char* argv[]) {
    me = argv[0];
    
    // parse input parameters
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-v" || arg == "--verbose") {
            __verbose = true;
        }
        else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) usage("input file");
            __input_name = argv[++i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) usage("output base name");
            __output_base = argv[++i];
        }
        else if (arg == "-r" || arg == "--radius") {
            if (i == argc-1) usage("radius");
            __radius = atof(argv[++i]);
        }
        else if (arg == "-d" || arg == "--density") {
            if (i == argc-1) usage("density");
            __density = atof(argv[++i]);
        }
        else if (arg == "-kn" || arg == "--knormal") {
            if (i >= argc-2) usage("spring constants");
            __k_load = atof(argv[++i]);
            __k_unload = atof(argv[++i]);
        }
        else if (arg == "-kr" || arg == "--kratio") {
            if (i == argc-1) usage("stiffness ratio");
            __k_ratio = atof(argv[++i]);
        }
        else if (arg == "-g" || arg == "--gamma") {
            if (i == argc-1) usage("gamma value");
            __gamma = atof(argv[++i]);
        }
        else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) usage("mu value");
            __mu = atof(argv[++i]);
        }
        else if (arg == "-f" || arg == "--frequency") {
            if (i == argc-1) usage("frequency");
            __f = atof(argv[++i]);
        }
        else if (arg == "-a" || arg == "--amplitude") {
            if (i == argc-1) usage("amplitude");
            __a = atof(argv[++i]);
        }
        else if (arg == "-dt" || arg == "--step") {
            if (i == argc-1) usage("integration step size");
            __dt = atof(argv[++i]);
        }
        else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) usage("tapping period");
            __T = atof(argv[++i]);
        }
        else if (arg == "-n" || arg == "--number") {
            if (i == argc-1) usage("number of taps");
            __ntaps = atoi(argv[++i]);
        }
        else if (arg == "-l" || arg == "--log") {
            if (i == argc-1) usage("log filename");
            __log_name = argv[++i];
        }
        else if (arg == "-of" || arg == "--outputf") {
            if (i == argc-1) usage("output frequency");
            __output_freq = atof(argv[++i]);
        }
        else usage("Invalid parameter: " + arg);
    }
    
    if (__input_name.empty()) usage("No input provided");
    if (__output_basename.empty()) usage("No output provided");
    
    std::ofstream __log_file;
    if (__verbose) {
        if (__log_name == "-")
            __log.rdbuf(std::cout.rdbuf());
        else {
            // dual log output to stdout and file
            __log_file.open(__log_name.c_str());
            if (__log_file.fail()) warning("Unable to open log file"); 
            
            __dual_sink_ptr = new dual_sink(std::cout, log_file);
            __dual_buffer_ptr = new dual_ostream_buffer(__dual_sink_ptr->rdbuf()));
            __log.rdbuf(*__dual_buffer_ptr);
        }
        
        // hijack std::cerr's ouput stream since all output must become part
        // of the requested log
        __cerr_buffer_ptr = std::cerr.rdbuf(__log.rdbuf());
    }
    else {
        __null_buffer_ptr = new null_ostream_buffer();
        __log.rdbuf(*__null_buffer_ptr);
        __cerr_buffer_ptr = 0;
    }

    __nthreads = 1;
#if _OPENMP
    __nthreads = omp_get_max_threads();
#endif
    
    // Load positions
    std::fstream __input_file(__input_name.c_str(), std::ios::in);
    if (__input_file.fail()) {
        warning("Unable to open " + __input_name, true);
    }
    boost::filesystem::path p(__input_name);
    std::string __filename = p.filename().string();
    std::string __extension = p.extension().string();
    if (__filename.substr(0, 3) == "in.") {
        import_data_liggghts(__input_file);
    }
    else if (__extension == ".txt") {
        import_data_plain(__input_file);
    }
    else {
        warning("Unrecognized input data format for " + __input_name, true);
    }
    __log << "File " << __filename << " imported successfully\n"
          << __particles.size() << " particles in input" << std::endl;
    const size_t nb_particles = __particles.size();
        
    // Set up traits
    __mass = 4./3.*M_PI*(__radius*__radius*__radius)*__density;
    particle_traits p_traits(__radius, __mass, __k_load, __k_unload,
                             __k_ratio, __gamma, __mu);
    boundary_traits b_traits();
    
    // Initialize cell grid
    nvis::timer __timer;
    nvis::uvec3 __res;
    __res[0] = floor(__bounds.size()[0] / __width);
    __res[1] = floor(__bounds.size()[1] / __width);
    __res[2] = std::min(res[0], res[1]);
    cell_grid __grid(__bounds, __periodic, __res, __particles);
    
    __log << "Cell grid with resolution " << __res 
          << " created in " << __timer.elapsed() << " s."
          << std::endl;
     
    // Initialize neighbor lists
    __neighbor_lists.resize(__particles.size());
    value_t rm = __radius + __skin;
    value_t rmsq = rm*rm;
    
    __timer.restart();
#pragma openmp parallel for
    for (size_t i=0 ; i<nb_particles ; ++i) {
        neighbor_list& nlist = __neighbor_lists[i];
        nlist = neighbor_list(i, __grid, __particles, rmsq);
    }
    __log << "Neighbor lists created in " << __timer.elapsed() << " s."
          << std::endl;
    
    __log << "Entering main loop" << std::endl;
    
    // 4. main loop
    //      4.a) compute per particle forces in parallel, update history
    //      4.b) add reciprocate forces in parallel
    //      4.c) take leapfrog step and update max displacement in parallel
    //      4.d) correct particle positions, if periodic boundaries
    //      4.e) compute diagnostic quantities, if needed
    //      4.f) update boundary position, if needed
    //      4.g) check max displacement against neighbor list update threshold
    //      4.h) update cell grid and neighbor list, if needed
    
    
    
    
    
    return 0;
}