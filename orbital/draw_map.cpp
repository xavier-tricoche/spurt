#include <iostream>
#include <fstream>
#include <vector>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>

#include "new_universal_poincare_map.hpp"
#include "cr3bp.hpp"
#include <misc/progress.hpp>

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

double K;
nvis::ivec2 res;
int maxi, nsamples;
std::string filename, input_name;
double eps;
double C;
double mu;
int seeding;
nvis::bbox2 bounds;
bool verbose = false;

typedef nvis::fixed_vector< double, 6 > vec6;


inline int pos(const nvis::vec2& x)
{
    const double& minx = bounds.min()[0];
    const double& miny = bounds.min()[1];
    const double& maxx = bounds.max()[0];
    const double& maxy = bounds.max()[1];

    int i = floor(res[0] * (x[0] - minx) / (maxx - minx));
    int j = floor(res[1] * (x[1] - miny) / (maxy - miny));
    if (i < 0 || i >= res[0] || j < 0 || j >= res[1]) return -1;
    return i + j*res[0];
}

inline void print_color(const nvis::vec3& col, int i, float* data)
{
    if (i < 0) return;
    data[3*i  ] += col[0];
    data[3*i+1] += col[1];
    data[3*i+2] += col[2];
}

void generate_seeds(std::vector< nvis::vec2 >& p)
{
    const double& minx = bounds.min()[0];
    const double& miny = bounds.min()[1];
    const double& maxx = bounds.max()[0];
    const double& maxy = bounds.max()[1];

    p.clear();
    switch (seeding) {
    case 0: {
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(minx + drand48()*(maxx - minx),
                                   miny + drand48()*(maxy - miny)));
        break;
    }
    case 1: {
        double dx = (maxx - minx) / (double)nsamples;
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(minx + (double)n*dx, 0.5*(miny + maxy)));
        break;
    }
    case 2: {
        double dy = (maxy - miny) / (double)nsamples;
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(0.5*(minx + maxx), miny + (double)n*dy));
        break;
    }
    case 3: {
        double dx = (maxx - minx) / (double)nsamples;
        double dy = (maxy - miny) / (double)nsamples;
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(minx + (double)n*dx, miny + (double)n*dy));
        break;
    }
    default: {
        std::cout << "unknown seeding pattern. using random sampling" << std::endl;
        seeding = 0;
        generate_seeds(p);
    }
    }
}

std::string me;
void displayUsageAndExit(const std::string& what) {
    if (what.size()) std::cerr << "ERROR: " << what << '\n';
    std::cerr
        << "USAGE: " << me << " [options]\n"
        << "DESCRIPTION: compute Poincare plot of circular restricted 3-body problem\n"
        << "OPTIONS:\n"
        << "   -h | --help                  Print this information\n"
        << "   -o | --output <string>       Output file name\n"
        << "   -r | --res <int> (x2)        Resolution\n"
        << "   -i | --iterations <int>      Number of iterations\n"
        << "   -b | --bounds <float> (x4)   Computation bounds (minx, miny, maxx, maxy)\n"
        << "   -e | --eps <float>           Integration precision\n"
        << "   -n | --samples <int>         Number of sample trajectories\n"
        << "      | --system <string>       System to consider\n"
        << "   -C <float>                   C constant\n"
        << "   -m | --mu <float>            mu constant\n"
        << "   -s | --seed <int>            Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal)\n"
        << "   -v | --verbose               Select verbose output\n";
    exit(1);
}

void readParams(const std::string& filename) {
    std::fstream input(filename, std::ios::in);
    while (!input.eof()) {
        std::string tmp, key;
        char c;
        input >> key;
        if (key == "mu") input >> mu;
        else if (key == "C") input >> C;
        else if (key == "xmin") input >> bounds.min()[0];
        else if (key == "ymin") input >> bounds.min()[1];
        else if (key == "xmax") input >> bounds.max()[0];
        else if (key == "ymax") input >> bounds.max()[1];
        else if (key == "eps") input >> eps;
        else if (key == "iter") input >> maxi;
        else if (key == "number") input >> nsamples;
        else if (key == "res") input >> res[0] >> res[1];
    }
}

int main(int argc, char* argv[])
{
    me = argv[0];
    bounds.min()[0] = 0.4;
    bounds.min()[1] = -2.5;
    bounds.max()[0] = 1.1;
    bounds.max()[1] = 2.5;
    C = 3.000;
    mu = Earth_Moon_mu;
    res = nvis::vec2(512, 512);
    eps = 1.0e-7;
    nsamples = 100;
    maxi = 1000;
    filename = "none";
    input_name = "none";

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") displayUsageAndExit("");
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) displayUsageAndExit("missing output name");
            filename = argv[++i];
        }
        else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) displayUsageAndExit("missing resolution");
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
        }
        else if (arg == "-C") {
            if (i == argc-1) displayUsageAndExit("missing C constant");
            C = atof(argv[++i]);
        }
        else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) displayUsageAndExit("missing mu constant");
            mu = atof(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) displayUsageAndExit("missing bounds");
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
        }
        else if (arg == "-i" || arg == "--iterations") {
            if (i == argc-1) displayUsageAndExit("missing number of iterations");
            maxi = atof(argv[++i]);
        }
        else if (arg == "-n" || arg == "--samples") {
            if (i == argc-1) displayUsageAndExit("missing number of samples");
            nsamples = atoi(argv[++i]);
        }
        else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) displayUsageAndExit("missing precision");
            eps = atof(argv[++i]);
        }
        else if (arg == "-s" || arg == "--seeding") {
            if (i == argc-1) displayUsageAndExit("missing precision");
            seeding = atoi(argv[++i]);
        }
        else if (arg == "-f" || arg == "--file") {
            if (i == argc-1) displayUsageAndExit("missing input filename");
            readParams(argv[++i]);
        }
        else if (arg == "--system") {
            if (i == argc-1) displayUsageAndExit("missing system name");
            std::string s(argv[++i]);
            if (s=="earth-moon" || s=="Earth-Moon") {
                mu = Earth_Moon_mu;
            }
            else if (s=="jupiter-europa" || s=="Jupiter-Europa") {
                mu = Jupiter_Europa_mu;
            }
        }
        else if (arg == "-v" || arg == "--verbose") {
            verbose=true;
        }
        else {
            displayUsageAndExit("unrecognized argument");
        }
    }

    if (nvis::all(bounds.min() == bounds.max()))
        displayUsageAndExit("Missing bounds");
    if (filename == "none")
        displayUsageAndExit("Missing output file");

    if (verbose) {
        std::cout << "parameters: resolution = " << res << '\n';
        std::cout << "bounds=\n" << bounds << std::endl;
    }

    if (nsamples == -1) nsamples = res[0];
    if (maxi == -1) maxi = res[0];

    float *hits = (float*)calloc(3 * res[0] * res[1], sizeof(float));

    srand48(987654321);

    std::vector< nvis::vec2 > seeds;
    generate_seeds(seeds);

    unsigned int count = 0;
    unsigned int pct, last_pct = 0;

    cr3bp dummy;
    planar_section<6>* section = new planar_section<6>(dummy.plane());
    int failed = 0;

    nvis::timer _timer;

    int nthreads=1;
#pragma omp parallel
    {
#ifdef _OPENMP
    nthreads=omp_get_num_threads();
#endif
    }
    if (nthreads==1) std::cout << "serial computation\n";
    else std::cout << "parallel computation with " << nthreads << " threads\n";

    spurt::ProgressDisplay progress;

    progress.start(nsamples);
    int counter=0;

#pragma omp parallel
    {
#ifdef _OPENMP
        const int thread_id = omp_get_thread_num();
#else
        const int thread_id = 0;
#endif
        cr3bp *field = new cr3bp(C, mu);
        new_universal_poincare_map<cr3bp, 6, planar_section<6> > _map(field, section);
        _map.precision(eps);

#pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < nsamples ; ++n) {
            ++count;
            nvis::vec2 x = seeds[n];
            nvis::vec3 color(drand48(), drand48(), drand48());
            print_color(color, pos(x), hits);
            vec6 y;
            try {
                y = field->unproject(x);
            }
            catch (...) {
                if (verbose) std::cout << "could not start at " << x << std::endl;
                ++failed;
                continue;
            }

#pragma omp atomic
            ++counter;

            std::vector< vec6 > __hits;
            try {
                _map.map(y, __hits, maxi);
            }
            catch (...) {
                std::cout << "caught exception at " << x << std::endl;
            }
            // std::cout << "everything went fine and we have " << __hits.size() << " hits" << std::endl;

            for (unsigned int i = 0 ; i < __hits.size() ; ++i) {
                print_color(color, pos(field->project(__hits[i])), hits);
            }

            pct = floor((double)count / (double)nsamples * 100.);
            if (thread_id==0) {
                progress.update(counter);
            }
//             if (thread_id == 0 && pct > last_pct)
//             if (false) {
//                 std::ostringstream os;
//                 os << count << " curves computed so far (" << pct << "%)\n";
//                 std::cout << os.str();
//                 last_pct = pct;
//             }
        }
    }
    progress.end();
    double delta_t = _timer.elapsed();
    std::cout << '\n';

    std::cout << failed << " integrations failed\n";
    std::cout << "integration of " << nsamples << " orbits took " << delta_t << " s. ("
         << static_cast<double>(nsamples)/delta_t << " Hz)\n";


    Nrrd *nout = nrrdNew();
    size_t size[3] = {3, (size_t)res[0], (size_t)res[1]};
    if (nrrdWrap_nva(nout, hits, nrrdTypeFloat, 3, size) ||
        nrrdSave(filename.c_str(), nout, NULL)) {
        std::cout << "ERROR while exporting FTLE file: " << biffGetDone(NRRD)
        << std::endl;
    }
    nrrdNuke(nout);

    return 0;
}
