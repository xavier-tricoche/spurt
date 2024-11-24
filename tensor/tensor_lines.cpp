#include <iostream>
#include <fstream>
#include <random>
#include <teem/nrrd.h>
#include <teem/unrrdu.h>

#include <math/types.hpp>
#include <math/inverse_transform.hpp>
#include <misc/progress.hpp>
#include <image/nrrd_field.hpp>
#include <flow/ftle.hpp>
#include <tensor/double_point.hpp>
#include <tensor/eigenvector_field.hpp>
#include <image/nrrd_wrapper.hpp>
#include <boost/random.hpp>

char* name_in;
char* name_out;
char* name_seed_density;
char* name_seed_points;
double length, sampling, dx, rel, min_length, max_length, d_length;
int dir, scaleLen, eigenid;
size_t nb_seeds;
char* tl_out;
int r;
float col[3];

typedef small_vector<double, 7> tensor_type;
typedef spurt::nrrd_field<tensor_type, 3, double> tensor_field_type;
typedef EigenvectorField<tensor_field_type> eigenvector_field_type;
typedef EigenvectorField<DoublePointLoad> dpl_eigenvectorfield_type;

void initialize(int argc, char* argv[], hestOpt* hopt)
{
    hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;

    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  0,  1,  &name_in,           "none",     "input file name, if none is provided, a synthetic double point load dataset will be procedurally generated");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,          NULL,       "output name");
    hestOptAdd(&hopt, "sd",     "seeding density",  airTypeString,  0,  1,  &name_seed_density, "none",     "seeding density volume name");
    hestOptAdd(&hopt, "sp",     "seed points",      airTypeString,  0,  1,  &name_seed_points,  "none",     "seed points file name");
    hestOptAdd(&hopt, "l",      "length",           airTypeDouble,  1,  1,  &length,            NULL,       "integration length for flow map computation");
    hestOptAdd(&hopt, "h",      "step size",        airTypeDouble,  1,  1,  &dx,                NULL,       "integration step size");
    hestOptAdd(&hopt, "n",      "nb seeds",         airTypeSize_t,  1,  1,  &nb_seeds,          NULL,       "number of seeds");
    hestOptAdd(&hopt, "random", NULL,               airTypeInt,     0,  0,  &r,                 NULL,       "random seeding");
    hestOptAdd(&hopt, "r",      "dpl rel distance", airTypeDouble,  0,  1,  &rel,               "0.5",      "relative distance between single point loads in procedural double point load model. Ignored if an input file is selected");
    hestOptAdd(&hopt, "e",      "eigenvector",      airTypeInt,     0,  1,  &eigenid,             "0",        "eigenvector field along which integration takes place");
    hestOptAdd(&hopt, "c",      "color",            airTypeFloat,   3,  3,  col,                "-1 -1 -1", "lines' color");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute tensor lines in eigenvector field of 3D symmetric second-order tensor field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct field_wrapper {
    field_wrapper(const DoublePointLoad& field)
        : dpl(new dpl_eigenvectorfield_type(field)),
          nrrd(0), procedural(true), bbox(field.bounds()) {}

    field_wrapper(const tensor_field_type& field)
        : dpl(0), nrrd(new eigenvector_field_type(field)),
          procedural(false), bbox(field.bounds()) {}

    spurt::vec3 interpolate(const spurt::vec3& x) const {
        if (procedural) {
            return dpl->evec(x, eigenid);
        } else {
            return nrrd->evec(x, eigenid);
        }
    }

    const spurt::bbox3& bounds() const {
        return bbox;
    }

    const dpl_eigenvectorfield_type*     dpl;
    const eigenvector_field_type*   nrrd;
    const bool                                   procedural;
    spurt::bbox3                                  bbox;
};

bool seeds_from_density(std::vector<spurt::vec3>& seeds, const std::string& name, size_t n)
{
    Nrrd* sampling_nrrd;
    spurt::bbox3 sampling_bounds;
    spurt::vec3 sampling_step;
    spurt::inverse_transform_sampling<float>* sampler;
    sampling_nrrd = spurt::nrrd_utils::readNrrd(name);
    int size = 1;
    if (sampling_nrrd->dim != 3) {
        std::cerr << "Provided NRRD seeding density file is not a 3D scalar volume. Ignored.\n";
        return false;
    }
    for (int i = 0 ; i < sampling_nrrd->dim ; ++i) {
        size *= sampling_nrrd->axis[i].size;
        if (std::isnan(sampling_nrrd->axis[i].min) || std::isinf(sampling_nrrd->axis[i].min)) {
            std::cerr << "Warning: no axis min provided for axis " << i << ". 0 assumed.\n";
            sampling_bounds.min()[i] = 0;
        } else {
            sampling_bounds.min()[i] = sampling_nrrd->axis[i].min;
        }
        if (std::isnan(sampling_nrrd->axis[i].spacing) || std::isinf(sampling_nrrd->axis[i].spacing)) {
            std::cerr << "Warning: no spacing information provided for axis " << i << ". 1 assumed.n";
            sampling_step[i] = 1.;
        } else {
            sampling_step[i] = sampling_nrrd->axis[i].spacing;
        }
    }
    std::vector<float> vals(size);
    spurt::nrrd_utils::to_vector(vals, sampling_nrrd);
    sampler = new spurt::inverse_transform_sampling<float>(vals);

    for (int count = 0 ; count < n ; ++count) {
        unsigned int id = sampler->sample();
        int i = id % sampling_nrrd->axis[0].size;
        int j = (id / sampling_nrrd->axis[0].size) % sampling_nrrd->axis[1].size;
        int k = id / (sampling_nrrd->axis[0].size * sampling_nrrd->axis[1].size);
        spurt::vec3 c(i, j, k);
        seeds.push_back(sampling_bounds.min() + (c * sampling_step));
    }

    return true;
}

bool seeds_from_file(std::vector<spurt::vec3>& seeds, const std::string& name, size_t n)
{
    seeds.clear();
    std::fstream in(name.c_str(), std::ios::in);
    if (in.fail()) {
        return false;
    }

    unsigned int N;
    in >> N;
    for (int i = 0 ; i < N && !in.eof() ; ++i) {
        spurt::vec3 x;
        in >> x[0] >> x[1] >> x[2];
        seeds.push_back(x);
    }

    if (seeds.size() > n) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(seeds.begin(), seeds.end(), g);
        seeds.resize(n);
    }

    return true;
}

void uniform_seeding(std::vector<spurt::vec3>& seeds, const spurt::bbox3& bounds, size_t n)
{
    seeds.clear();
    spurt::vec3 diameter = bounds.size();
    double ref = *std::min_element(&diameter[0], &diameter[3]);
    int fac = 1;
    for (int i = 0 ; i < 3 ; ++i) {
        fac *= diameter[i] / ref;
    }
    float refn = pow((float)n / fac, 1. / 3.);
    int nsamples[3];
    spurt::vec3 step;
    for (int i = 0 ; i < 3 ; ++i) {
        nsamples[i] = (int)floor(refn * diameter[i] / ref);
        step[i] = std::max(1., diameter[i] / (float)(nsamples[i] - 1));
    }
    int npoints = nsamples[0] * nsamples[1] * nsamples[2];

    for (int l = 0 ; l < npoints ; ++l) {
        int i = l % nsamples[0];
        int j = (l / nsamples[0]) % nsamples[1];
        int k = l / (nsamples[0] * nsamples[1]);
        spurt::vec3 c(i, j, k);
        spurt::vec3 x = bounds.min() + (c * step);
        seeds.push_back(x);
    }
}

void random_seeding(std::vector<spurt::vec3>& seeds, const spurt::bbox3& bounds, size_t n)
{
    srand48(time(0));
    seeds.clear();
    spurt::vec3 diameter = bounds.size();
    for (int i = 0 ; i < n ; ++i) {
        spurt::vec3 c(drand48()*diameter[0], drand48()*diameter[1], drand48()*diameter[2]);
        spurt::vec3 x = bounds.min() + c;
        seeds.push_back(x);
    }
}

int main(int argc, char* argv[])
{
    hestOpt* hopt;
    initialize(argc, argv, hopt);

    if (eigenid < 0 || eigenid > 2) {
        std::cerr << "ERROR: invalid eigenvector field selected" << std::endl;
        hestUsage(stderr, hopt, argv[0], 0);
        return -1;
    }

    bool procedural = !strcmp(name_in, "none");
    field_wrapper*   efield         = 0;
    DoublePointLoad* dpl            = 0;
    tensor_field_type*   nrrd_tensor    = 0;
    if (procedural) {
        dpl = new DoublePointLoad(50, 50, 20, rel);
        dpl->set_check_inside(false);
        efield = new field_wrapper(*dpl);
        std::cerr << "generating a synthetic double point load tensor field" << std::endl;
    } else {
        nrrd_tensor = new tensor_field_type(name_in);
        efield = new field_wrapper(*nrrd_tensor);
        std::cerr << "processing nrrd file: " << name_in << std::endl;
    }

    std::vector<spurt::vec3> seeds;
    bool pdf_seeding = strcmp(name_seed_density, "none");
    bool file_seeding = strcmp(name_seed_points, "none");

    if (file_seeding) {
        std::cerr << "reading seed points from file " << name_seed_points << "\n";
        seeds_from_file(seeds, name_seed_points, nb_seeds);
    } else if (pdf_seeding) {
        std::cerr << "choosing seed points based on density file " << name_seed_density << '\n';
        seeds_from_density(seeds, name_seed_density, nb_seeds);
    } else if (r) {
        std::cerr << "choosing seed points at random\n";
        random_seeding(seeds, efield->bounds(), nb_seeds);
    } else {
        std::cerr << "uniform seeding\n";
        uniform_seeding(seeds, efield->bounds(), nb_seeds);
    }

    int lastpct = -1;

    std::vector<std::vector<spurt::vec3> > curves;
    int nverts = 0;
    int nprocessed = 0;
    spurt::timer timer;
    for (int n = 0 ; n < seeds.size() ; ++n) {
        int pct = 100 * n / seeds.size();
        if (pct > lastpct) {
            lastpct = pct;
            std::cerr << '\r' << pct << "% completed in "
                      << timer.elapsed()
                      << "s.                 \r"
                      << std::flush;
        }

        spurt::vec3 seed = seeds[n];

        for (int dir = 0 ; dir < 2 ; ++dir) {
            int error = -3;
            double h = (dir ? dx : -dx);
            curves.push_back(std::vector<spurt::vec3>());
            std::vector<spurt::vec3>& steps = curves.back();
            try {
                double length_io = length;
                spurt::vec3 z = ftle::eigen_flow_map(*efield, seed, h, length_io, error, steps);
                if (spurt::any(spurt::isinvalid(z))) {
                    continue;
                }
            } catch (...) {
                continue;
            }
            nverts += curves.back().size();
        }
    }

    std::cout << "\ntotal computation time for tensor lines was " << timer.elapsed() << '\n';

    std::fstream out(name_out, std::ios::out);
    out << "# vtk DataFile Version 2.0\n"
        << "Tensor lines computed in " << name_in << '\n'
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << nverts << " float\n";
    int nlines = 0;
    int nptsinlines = 0;
    for (int i = 0 ; i < curves.size() ; ++i) {
        for (int j = 0 ; j < curves[i].size() ; ++j) {
            out << curves[i][j][0] << " " << curves[i][j][1] << " " << curves[i][j][2] << '\n';
        }
        if (curves[i].size() > 2) {
            ++nlines;
            nptsinlines += curves[i].size();
        }
    }
    out << "LINES " << nlines << " " << nptsinlines + nlines << '\n';
    int offset = 0;
    for (int i = 0 ; i < curves.size() ; ++i) {
        if (curves[i].size() < 3) {
            offset += curves[i].size();
            continue;
        }
        out << curves[i].size();
        for (int j = 0 ; j < curves[i].size() ; ++j) {
            out << " " << offset + j;
        }
        out << '\n';
        offset += curves[i].size();
    }
    if (col[0] > 0) {
        out << "POINT_DATA " << nverts << '\n';
        out << "COLOR_SCALARS " << nverts << " 3\n";
        for (int i=0 ; i<nverts ; ++i) {
            out << col[0] << " " << col[1] << " " << col[2] << '\n';
        }
    }
    out.close();

    if (procedural) {
        delete dpl;
    } else {
        delete nrrd_tensor;
    }

    return 0;
}
