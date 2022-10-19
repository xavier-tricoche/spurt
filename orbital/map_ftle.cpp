#include <flow/ftle.hpp>
#include <poincare/metric.hpp>

#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <boost/format.hpp>
#include <boost/limits.hpp>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <util/wall_timer.hpp>

#include <math/math.hpp>
#include "new_universal_poincare_map.hpp"
#include "cr3bp.hpp"
#include <format/png_io.hpp>
#include <format/filename.hpp>
#include <misc/meta_utils.hpp>
#include <misc/misc_helper.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>

// #include <data/raster.hpp>

#include <teem/nrrd.h>

#ifndef NO_TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>
#endif

#ifndef NO_TBB
tbb::mutex output_mutex;
tbb::mutex progress_mutex;
tbb::mutex orbit_mutex;
tbb::atomic<size_t> progress_counter;
#else
size_t progress_counter;
#endif

void write_to_ostream(std::ostream& os, const std::string& str) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(output_mutex);
    #endif
        os << str << '\n';
    }
}

void update_progress(xavier::ProgressDisplay& progress) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(progress_mutex);
    #endif
        progress.update(progress_counter);
    }
}

double minx, maxx, miny, maxy;

namespace nvis {
    typedef nvis::fixed_vector<double, 6> vec6;
}

using namespace nvis;

double eps;
int nx, ny, maxit, deltait;
double hx, hy;
std::vector < std::vector < bool > > valid;
char *outs, *sys;
bool is_periodic = false;
double C;
double mu;
nvis::bounding_box<nvis::vec2> bounds;
bool extract_tangent_curves = false;

typedef nvis::fixed_vector<double, 6>    vec6;
typedef xavier::raster_grid<2, double> grid_t;
typedef grid_t::coord_type coord_t;

template<typename Data_>
using raster_t=xavier::raster_data<Data_, 2>;

inline double distance(const nvis::vec2& p0, const nvis::vec2& p1)
{
    return nvis::norm(p1 - p0);
}

void initialize(int argc, char* argv[])
{
    hestOpt *hopt = NULL;
    hestParm *hparm;
    airArray *mop;
    char *me;

    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "o",    "output",          airTypeString, 1, 1, &outs,  NULL,      "output name");
    hestOptAdd(&hopt, "nx",   "resolution x",    airTypeInt,    0, 1, &nx,    "512",  "sampling resolution along x");
    hestOptAdd(&hopt, "ny",   "resolution y",    airTypeInt,    0, 1, &ny,    "512",  "sampling resolution along y");
    hestOptAdd(&hopt, "maxi", "max iterations",  airTypeInt,    1, 1, &maxit, NULL,      "max number of map iterations");
    hestOptAdd(&hopt, "stride", "iteration stride",  airTypeInt,    1, 1, &deltait, NULL,      "Output iteration stride");
    hestOptAdd(&hopt, "e",    "epsilon",         airTypeDouble, 0, 1, &eps,   "1.0e-6",  "integration accuracy");
    hestOptAdd(&hopt, "minx", "min x coord",     airTypeDouble, 0, 1, &minx,  "-1",      "min x coord. of bounding box");
    hestOptAdd(&hopt, "miny", "min y coord",     airTypeDouble, 0, 1, &miny,  "-1.3335", "min y coord. of bounding box");
    hestOptAdd(&hopt, "maxx", "max x coord",     airTypeDouble, 0, 1, &maxx,  "0",       "max x coord. of bounding box");
    hestOptAdd(&hopt, "maxy", "max y coord",     airTypeDouble, 0, 1, &maxy,  "1.3355",  "max y coord. of bounding box");
    hestOptAdd(&hopt, "C",    "Jacobi constant", airTypeDouble, 0, 1, &C,     "4.5",     "Jacobi constant");
    hestOptAdd(&hopt, "mu",   "mu constant",     airTypeDouble, 0, 1, &mu,    "0.5",     "mu constant");
    hestOptAdd(&hopt, "s",    "system",          airTypeString, 0, 1, &sys,   "earth",   "system name");
    hestOptAdd(&hopt, "t",    "tangent",         airTypeBool,   0, 1, &extract_tangent_curves, "0", "Extract tangent curves");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE value of circular restricted 3-body problem",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline int pos(const nvis::vec2& x)
{
    int i = floor((x[0] - minx) / hx);
    int j = floor((x[1] - miny) / hy);
    return i + j*nx;
}

inline nvis::vec2 pos(int i, int j)
{
    return nvis::vec2(minx + (double)i*hx, miny + (double)j*hy);
}

// Looks for zero crossings of ydot along the edges
// of a 2d cell.
inline bool is_tangent(const nvis::vec4& dys) {
    static const int ids[5] = { 0, 1, 2, 3, 0 };
    for (int i=0; i<4; ++i) {
        if (dys[i]*dys[i+1]<0) {
            std::cout << "tangent!!\n";
            return true;
        }
    }
    return false;
}

double lmax(size_t n, const raster_t<vec6>& pos, const raster_t<bool>& _valid, const cr3bp& field)
{
    if (!_valid[n]) return -1; // invalid value

    size_t i = n % nx;
    size_t j = n / nx;
    std::vector<size_t> indices[2];
    if (i==0) {
        indices[0].resize(2);
        indices[0][0] = n;
        indices[0][1] = n+1;
    }
    else if (i==nx-1) {
        indices[0].resize(2);
        indices[0][0] = n-1;
        indices[0][1] = n;
    }
    else {
        indices[0].resize(3);
        indices[0][0] = n-1;
        indices[0][1] = n;
        indices[0][2] = n+1;
    }
    if (j==0) {
        indices[1].resize(2);
        indices[1][0] = n;
        indices[1][1] = n+nx;
    }
    else if (j==ny-1) {
        indices[1].resize(2);
        indices[1][0] = n-nx;
        indices[1][1] = n;
    }
    else {
        indices[1].resize(3);
        indices[1][0] = n-nx;
        indices[1][1] = n;
        indices[1][2] = n+nx;
    }
    // if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || !_valid[n]) return -1.0;

    static double dh[2] = { hx, hy };
    nvis::vec2 dxi[2];
    // int indices[2][3] = { { n - 1, n, n + 1 }, { n - nx, n, n + nx } };
    for (unsigned int dim = 0 ; dim < 2 ; ++dim) {
        int min = -1;
        int max = min;
        for (int i = 0 ; i < indices[dim].size() ; ++i) {
            if (_valid[indices[dim][i]]) {
                if (min < 0) min = i; // this is the first valid index
                else max = i;
            }
        }
        if (max <= min) return -1.0;
        dxi[dim] = 1. / (dh[dim] * (double)(max - min)) *
            distance(field.project(pos[indices[dim][min]]), field.project(pos[indices[dim][max]]));
    }
    double a = nvis::inner(dxi[0], dxi[0]);
    double b = nvis::inner(dxi[0], dxi[1]);
    double c = nvis::inner(dxi[1], dxi[1]);

    double lmaj = 0.5 * (a * a + 2 * b * b + c * c + (a + c) * sqrt((a - c) * (a - c) + 4 * b * b));

    return lmaj;
}

template<typename T>
void output(const std::string& name, const raster_t<T>& data) {
    typedef typename xavier::data_traits<T>::value_type value_t;

    std::cout << "exporting " << name << std::endl;
    std::ofstream outf(name.c_str(), std::ofstream::binary);
    std::ostringstream os;
    std::string sizestr, spacingstr, minstr;
    size_t sz=xavier::data_traits<T>::size();
    if (sz>1) os << sz << " ";
    os << data.grid().resolution()[0] << " "
       << data.grid().resolution()[1];
    sizestr=os.str();
    os.str("");
    if (sz>1) os << "nan ";
    os << data.grid().spacing()[0] << " "
       << data.grid().spacing()[1];
    spacingstr=os.str();
    os.str("");
    if (sz>1) os << "nan ";
    os << data.grid().bounds().min()[0] << " "
       << data.grid().bounds().min()[1];
    minstr = os.str();
    os.str("");
    os << "NRRD0001\n"
        << "# Complete NRRD file format specification at:\n"
        << "# http://teem.sourceforge.net/nrrd/format.html\n"
        << "type: " << xavier::type2string<value_t>::type_name() << "\n"
        << "dimension: " << (sz>1 ? 3 : 2) << "\n"
        << "sizes: " << sizestr << '\n'
        << "spacings: " << spacingstr << '\n'
        << "axis mins: " << minstr << '\n'
        << "endian: little\n"
        << "encoding: raw\n\n";
    outf.write(os.str().c_str(), os.str().size());
    size_t size=sz*data.size()*sizeof(T);
    outf.write((char*)&data[0], size);
    outf.close();
    std::cout << name << " has been exported\n";
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);

    if (sys != NULL) {
        std::string name(sys);
        name = xavier::lower_case(name);
        if (name=="earth" || name=="earth_moon" || name=="earth-moon")
            mu=Earth_Moon_mu;
        else if (name=="jupiter" || name=="jupiter_europa" ||
                 name=="jupiter-europa")
            mu=Jupiter_Europa_mu;
        else {
            std::cerr << "Warning: System name " << name << " unrecognized\n";
        }
    }

    std::cout << "max iterations = " << maxit
              << ", eps = " << eps << '\n';
    hx=(maxx-minx)/(nx-1);
    hy=(maxy-miny)/(ny-1);
    std::cout << "sampling resolution is: hx = " << hx << ", hy = " << hy << std::endl;
    std::cout << "mesh resolution is " << nx << " x " << ny << std::endl;

    bounds.min() = nvis::vec2(minx, miny);
    bounds.max() = nvis::vec2(maxx, maxy);

    cr3bp *field = new cr3bp(C, mu);

    valid.resize(maxit);
    for (unsigned int p = 0 ; p < maxit ; ++p) {
        valid[p].resize(nx*ny);
        std::fill(valid[p].begin(), valid[p].end(), false);
    }

    planar_section<6>* plane = new planar_section<6>(field->plane());
    grid_t grid(nvis::ivec2(nx, ny), bounds);
    raster_t<nvis::fvec2> ftle(grid);
    raster_t<nvis::vec6> pos_fwd(grid), pos_bwd(grid);
    raster_t<bool> valid_fwd(grid, true), valid_bwd(grid, true);
    raster_t<double> times_fwd(grid, 0.), times_bwd(grid, 0.);
    raster_t<unsigned char> tangent(grid, 0);
    raster_t<float> dycomp(grid, 0);

    // initialize positions and determine tangency locations
    const nvis::vec2& spc=grid.spacing();

    // vertices of dual grid cell
    nvis::vec6 dx(0), dxdot(0);
    dx[0]=spc[0]/2;
    dxdot[3]=spc[1]/2;
    std::vector<nvis::vec6> dp(4);
    dp[0]=-1*dx-dxdot;
    dp[1]=dx-dxdot;
    dp[2]=-1*dp[0];
    dp[3]=-1*dp[1];

    xavier::ProgressDisplay progress(true);
    progress.fraction_on();
    progress_counter = 0;
    progress.begin(nx*ny, "Checking initial conditions");
#ifdef NO_TBB
    for (int n = 0 ; n < nx*ny ; ++n) { // n is a seed point index
#else
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nx*ny),
                  [&](tbb::blocked_range<size_t> r) {
    for (size_t n=r.begin(); n!=r.end(); ++n) {
#endif
        ++progress_counter;
        update_progress(progress);
        try {
            pos_bwd[n] = pos_fwd[n] = field->unproject(grid(grid.coordinates(n)));
        }
        catch (...) {
            valid_fwd[n] = valid_bwd[n] = false;
            tangent[n]=0;
        }
        if (valid_fwd[n]) {
            try {
                nvis::vec6 f=(*field)(0, pos_fwd[n]);
                dycomp[n]=f[1];
            }
            catch(...) {}
            bool valid=true;
            // check for tangency in dual cell centered at current position
            std::vector<nvis::vec6> p(4, pos_fwd[n]);
            nvis::vec4 dy(4);
            for (int i=0; i<4; ++i) {
                p[i]+=dp[i];
                nvis::vec6 f;
                try {
                    f=(*field)(0, p[i]);
                }
                catch(...) {
                    valid=false;
                    break;
                }
                dy[i]=f[1]; // dX/dy
            }
            if (valid) tangent[n]=is_tangent(dy) ? 2 : 1;
        }
    }
#ifndef NO_TBB
    });
#endif
    progress.end();
    {
        std::ostringstream os;
        os << outs << "-C=" << C
           << "-bounds=[" << minx << "," << maxx
           << "]x[" << miny << "," << maxy << "]-tangent.nrrd";
        output(os.str(), tangent);
        os.str("");
        os << outs << "-C=" << C
           << "-bounds=[" << minx << "," << maxx
           << "]x[" << miny << "," << maxy << "]-normal.nrrd";
        output(os.str(), dycomp);
    }
    int last_nlines = 0;

    std::vector<int> random_points(100);
    for (int i=0; i<100; ++i) {
        random_points[i] = rand() % (nx*ny);
    }

    for (unsigned int p = 0 ; p < maxit ; ++p) {

        std::cerr << "\niteration " << p << '\n' << std::flush;

        std::fill(ftle.begin(), ftle.end(), nvis::fvec2(0,0));

        {
            new_universal_poincare_map<cr3bp, 6, planar_section<6> > my_map(field, plane);
            my_map.precision(eps);

            progress.start(nx*ny, "Integrating orbits");
            progress_counter = 0;
#ifdef NO_TBB
            for (int n = 0 ; n < nx*ny ; ++n) { // n is a seed point index
#else
            tbb::parallel_for(tbb::blocked_range<size_t>(0, nx*ny),
                  [&](tbb::blocked_range<size_t> r) {
            for (size_t n=r.begin(); n!=r.end(); ++n) {
#endif
                ++progress_counter;
                update_progress(progress);

                std::vector< vec6 > hits;
                std::vector< double > times;
                if (valid_fwd[n]) {
                    vec6 x = pos_fwd[n];
                    try {
                        hits.clear();
                        times.clear();
                        my_map.map(x, hits, times, 1);
                        pos_fwd[n] = hits[0];
                        times_fwd[n] += times[0];
                    }
                    catch (...) {
                        valid_fwd[n] = false;
                    }
                }

                if (valid_bwd[n]) {
                    vec6 x = pos_bwd[n];
                    try {
                        hits.clear();
                        times.clear();
                        my_map.map(x, hits, times, -1);
                        pos_bwd[n] = hits[0];
                        times_bwd[n] += times[0];
                    }
                    catch (...) {
                        valid_bwd[n] = false;
                    }
                }
            }
#ifndef NO_TBB
            });
#endif
        }
        progress.end();

        progress.start(nx*ny, "Computing FTLE");
        progress_counter = 0;
        {
#ifdef NO_TBB
            for (int n = 0 ; n < nx*ny ; ++n) { // n is a seed point index
#else
            tbb::parallel_for(tbb::blocked_range<size_t>(0, nx*ny),
                              [&](tbb::blocked_range<size_t> r) {
            for (size_t n=r.begin(); n!=r.end(); ++n) {
#endif
                ++progress_counter;
                update_progress(progress);

                if (!valid_fwd[n]) continue;
                if (tangent[n]<2) {
                    double ftle_fwd = log(std::max(lmax(n, pos_fwd, valid_fwd, *field), 1.)) / times_fwd[n];
                    double ftle_bwd = log(std::max(lmax(n, pos_bwd, valid_bwd, *field), 1.)) / fabs(times_bwd[n]);
                    ftle[n]=nvis::fvec2(ftle_fwd, ftle_bwd);
                }
            }
#ifndef NO_TBB
            });
#endif
        }
        std::cout << progress << '\n';

        std::ostringstream os;
        os << xavier::filename::remove_extension(outs)
           << "-ftle-C=" << C
           << "-bounds=[" << minx << "," << maxx
           << "]x[" << miny << "," << maxy << "]"
           << "-eps=" << eps
           << "-step_" << p+1 << "_of_" << maxit << ".nrrd";
        output(os.str(), ftle);
    }

    return 0;
}
