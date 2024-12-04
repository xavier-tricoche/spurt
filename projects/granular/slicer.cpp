#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <teem/nrrd.h>
#include <assert.h>
#include <vector>
#include <sstream>
#include <iomanip>
#include <image/nrrd_wrapper.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif
// spurt
#include <image/nrrd_wrapper.hpp>

inline double cubed(double x)
{
    return x*x*x;
}

// compute volume fraction of a sphere of radius r at height y that
// falls in slab [lo, hi]
inline double fraction(double r, double y, double lo, double hi)
{
    // upper and lower bound in local coordinates  of sphere
    double y0 = std::max(lo-y, -r);
    double y1 = std::min(hi-y, r);
    if (y1 <= y0) {
        return 0.;
    }
    return M_PI*(r*r*(y1 - y0) + (cubed(y1) - cubed(y0))/3.);
}

// compute the volume fraction of a particle of radius r at height y
// that falls in each slab of thickness dy. y0 is the floor height.
inline void slice(double r, double y, double y0, double dy,
                    std::vector<std::pair<int, double> >& fracs)
{
    // local coordinates
    double z = y-y0;
    // local bounds
    double zmin = z-r;
    double zmax = z+r;
    // range of intersected slabs
    int imin = (int)(zmin/dy);
    int imax = (int)(zmax/dy);
    
    typedef std::pair<int, double> pair_type;
    for (int i=imin ; i<=imax ; ++i) {
        fracs.push_back(pair_type(i, fraction(r, y, dy*i, dy*(i+1))));
    }
}

std::string me, in_name, out_name;
double dy=0.01, minx=0, miny=0, minz=0, maxx=0.24, maxz=0.24, maxy=1.0, radius=0.01;
int nb=3456;
bool periodic = true;

void printUsageAndExit(const std::string& offending="",
                       bool doExit = true )
{
    if (offending != "") {
        std::cerr << "ERROR: " << offending << std::endl;
    }
    std::cerr
            << "Usage  : " << me     << " [parameters] [options]\n"
            << "Parameters:\n"
            << "    -i  | --input <path>              Input file name\n"
            << "    -o  | --output <basename>         Output file name\n"
            << "Options:\n"
            << "    -h  | --help                      Print this information\n"
            << "    -r  | --radius <float>            Particle radius\n"
            << "    -s  | --step <float>              Slicing step size\n"
            << "    -n  | --number <int>              Number of particles in system\n"
            << "    -bx | --boundsx <float> <float>   Domain bounding box along x-axis\n"
            << "    -by | --boundsy <float> <float>   Domain bounding box along y-axis\n"
            << "    -bz | --boundsz <float> <float>   Domain bounding box along z-axis\n"
            << "    -p  | --periodic <bool>           Domain is periodic in X/Z\n"
            << std::endl;
            
    if (doExit) {
        exit(1);
    }
}

std::vector<float> data;
std::fstream input;
bool isnrrd;

// import a single line of data from a NJIT file
inline bool read_line_coord(float& t, float& x, float& y, float& z,
                            std::fstream& file)
{
    int i;
    file >> i >> x >> y >> z >> t;
    return !file.eof();
}

// import the height of each particle at a given time step
inline bool getY(float& y, int ts, int i)
{
    if (isnrrd) {
        y=data[3*nb*ts+3*i+1];
    } else {
        float t, x, z;
        if (!read_line_coord(t, x, y, z, input)) {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    me = argv[0];
    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                printUsageAndExit();
            }
            in_name = argv[++i];
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit();
            }
            out_name = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            printUsageAndExit();
        } else if (arg == "-s" || arg == "--step") {
            if (i == argc-1) {
                printUsageAndExit();
            }
            dy = atof(argv[++i]);
        } else if (arg == "-n" || arg == "--number") {
            if (i == argc-1) {
                printUsageAndExit();
            }
            nb = atoi(argv[++i]);
        } else if (arg == "-bx" || arg == "--boundx") {
            if (i >= argc-2) {
                printUsageAndExit();
            }
            minx = atof(argv[++i]);
            maxx = atof(argv[++i]);
        } else if (arg == "-by" || arg == "--boundy") {
            if (i >= argc-2) {
                printUsageAndExit();
            }
            miny = atof(argv[++i]);
            maxy = atof(argv[++i]);
        } else if (arg == "-bz" || arg == "--boundz") {
            if (i >= argc-2) {
                printUsageAndExit();
            }
            minz = atof(argv[++i]);
            maxz = atof(argv[++i]);
        } else if (arg == "-p" || arg == "--periodic") {
            if (i == argc-1) {
                printUsageAndExit();
            }
            periodic = atoi(argv[++i]);
        } else if (arg == "-r" || arg == "--radius") {
            if (i == argc-1) {
                printUsageAndExit();
            }
            radius = atof(argv[++i]);
        }
    }
    
    size_t found = in_name.find_last_of('.');
    std::string ext = in_name.substr(found+1);
    int nbdim=3;
    if (ext == "nrrd") {
        Nrrd* nin = nrrdNew();
        if (nrrdLoad(nin, in_name.c_str(), NULL)) {
            printUsageAndExit(biffGetDone(NRRD));
        }
        spurt::to_vector<float>(data, nin);
        isnrrd = true;
        nbdim = nin->dim;
        nrrdNuke(nin);
    } else if (ext == "txt") {
        input.open(in_name.c_str(), std::ios::in);
        if (!input) {
            std::cerr << me << ": unable to open " << in_name;
            printUsageAndExit();
        }
        isnrrd = false;
    } else {
        printUsageAndExit("unrecognized input file extension");
    }
    
    size_t nbslices = (maxy - miny)/dy;
    float slicevol = (maxx-minx)*dy*(maxz-minz);
    std::vector<float> particles(nb);
    for (int counter=0 ; true ; ++counter) {
        int i;
        for (i=0 ; i<nb ; ++i) {
            if (!getY(particles[i], counter, i)) {
                break;
            }
        }
        if (i < nb) {
            std::cerr << "time step prematurily interrupted. exit" << std::endl;
            exit(-1);
        }
        
        std::vector<std::vector<std::pair<int, double> > > fractions(nb);
        
#pragma openmp parallel for
        for (int n=0 ; n<nb ; ++n) {
            slice(radius, particles[n], miny, dy, fractions[n]);
        }
        
        float* data = (float*)calloc(2*nbslices, sizeof(float));
        for (int n=0 ; n<nb ; ++n) {
            const std::vector<std::pair<int, double> >& frac = fractions[n];
            for (int m=0 ; m<frac.size() ; ++m) {
                data[2*frac[m].first+1] += frac[m].second;
            }
        }
        float measured_vol = 0;
        for (int m=0 ; m<nbslices ; ++m) {
            data[2*m] = m*dy;
            measured_vol += data[2*m+1];
            data[2*m+1] /= slicevol;
        }
        float actual_vol = nb * 4./3. * M_PI * radius*radius*radius;
        
        std::ostringstream os;
        if (nbdim == 3) {
            os << out_name << "_" << std::setw(6) << std::setfill('0') << counter << ".nrrd";
        } else {
            os << out_name;
        }
        Nrrd* nout = nrrdNew();
        size_t size[] = {2, nbslices};
        if (nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, size) ||
                nrrdSave(os.str().c_str(), nout, NULL)) {
            std::cerr << biffGetDone(NRRD) << std::endl;
            exit(-1);
        }
        std::cerr
                << "exported:\t" << os.str()
                << " (" << measured_vol/actual_vol*100 << "%)" << std::endl;
        nrrdNuke(nout);
        if (nbdim == 2) {
            break;
        }
    }
    
    return 0;
}

