#ifndef __FTLE_HPP__
#define __FTLE_HPP__

#include <math/types.hpp>
#include <misc/progress.hpp>
#include <math/types.hpp>
#include <sstream>
#include <string>
#include <misc/misc_helper.hpp>

struct file_name {

    file_name(const std::string& base) : basename(base) {}

    std::string make(const std::string& what, int id, float length = 0) {
        os.clear();
        os.str("");
        os << basename << "-" << what << "--" << id;
        if (length > 0) {
            os << "--l=" << length;
        }
        os << ".nrrd";
        return os.str();
    }

    std::ostringstream os;
    std::string basename;
};


namespace {

template< typename EigenvectorField >
struct OrientationChecker {
    OrientationChecker(const EigenvectorField& field, const spurt::vec3& x)
        : _field(field) {
        _ref_ev = _field.interpolate(x);
    }
    bool consistent(const spurt::vec3& x) const {
        spurt::vec3 ev = _field.interpolate(x);
        return (spurt::inner(_ref_ev, ev) > 0);
    }

    const EigenvectorField& _field;
    spurt::vec3 _ref_ev;
};

// helper function to compute flow map derivative on bi-directional flow
template < typename FlowMap, typename Checker >
void central_diff(spurt::vec3 ans[2], const FlowMap flowmaps[2],
                  const Checker& checker, const spurt::ivec3& coord, int dim)
{
    int i = coord[dim];
    int imax = flowmaps[0].grid().resolution()[dim];
    double h = flowmaps[0].grid().spacing()[dim];

    spurt::ivec3 ids[2] = { coord, coord }; // "left" and "right" neighbors
    if (i > 0) {
        ids[0][dim] = i - 1;
    }
    if (i < imax - 1) {
        ids[1][dim] = i + 1;
    }
    assert(ids[0][dim] < ids[1][dim]);

    double dx = (double)(ids[1][dim] - ids[0][dim]) * h;

    int fld_id[2] = {0, 0};
    for (int i = 0 ; i < 2 ; ++i) {
        if (!checker.consistent(flowmaps[0].grid()(ids[i]))) {
            fld_id[i] = 1;
        }
    }

    // forward derivative
    ans[0] = flowmaps[fld_id[1]](ids[1]) - flowmaps[fld_id[0]](ids[0]);
    ans[1] = flowmaps[1-fld_id[1]](ids[1]) - flowmaps[1-fld_id[0]](ids[0]);
    ans[0] /= dx;
    ans[1] /= dx;
}
}

namespace ftle {

template< typename Field >
spurt::vec3 flow_map(const Field& field, const spurt::vec3& seed, double dx, double length)
{
    spurt::vec3 x = seed, lastx = seed;
    double l = 0, lastl = 0;
    while (l < length) {
        try {
            spurt::vec3 step = dx * field.interpolate(x);
            lastx = x;
            lastl = l;
            x += step;
            if (spurt::norm(step) < 1.0e-9) {
                return x;
            }
            l += spurt::norm(step);
        } catch (...) {
            return lastx;
        }
    }
    double u = (length - lastl) / (l - lastl);
    return (1. - u)*lastx + u*x;
}

template< typename Field >
spurt::vec3 eigen_flow_map(const Field& field, const spurt::vec3& seed, double dt, double& delta_T, int& error,
                          spurt::vec3& start_dir, std::vector<spurt::vec3>& steps, int nmax=10000)
{
    spurt::vec3 x = seed, lastx = seed;
    double t = 0, lastt = 0;

    steps.clear();
    steps.push_back(seed);

    try {
        // initialize integration direction
        spurt::vec3 ref_dir = field.interpolate(x);
        if (spurt::norm(start_dir) && spurt::inner(ref_dir, start_dir) < 0) {
            dt *= -1;
            ref_dir *= -1;
        }

        ref_dir /= spurt::norm(ref_dir); // ref_dir is always normalized

        int n=0;
        while (fabs(t) < delta_T && n<nmax) {
            spurt::vec3 step = dt * field.interpolate(x);
            // check orientation consistency
            double dot = spurt::inner(step, ref_dir);
            step *= spurt::sign(dot);
            start_dir = step;

            // std::cerr << "at " << x << ", t = " << t << ", step = " << step << std::endl;

            double step_length = spurt::norm(step);
            if (dot / step_length < 0.5) {
                // halve step size in case of sharp turn
                step *= 0.5;
                step_length *= 0.5;
            }
            if (step_length < 1.0e-9) {
                start_dir /= spurt::norm(start_dir);
                return x;
            }
            // update state machine
            ref_dir = step / step_length;
            lastx = x;
            lastt = t;
            // advance
            x += step;
            t += dt;
            steps.push_back(x);
        }
        if (t == 0) {
            delta_T = 0;
            throw std::runtime_error("1000"); // we went nowhere
        }
        start_dir /= spurt::norm(start_dir);
        error = 0; // everything went fine
        double u = (delta_T - lastt) / (t - lastt);
        steps.push_back((1. - u)*lastx + u*x);
        return steps.back();
    } catch (std::runtime_error& e) {
        error = atoi(e.what());
        delta_T = lastt;
        start_dir /= spurt::norm(start_dir);
        return lastx;
    }
}

template< typename Field >
spurt::vec3 eigen_flow_map(const Field& field, const spurt::vec3& seed, double dx, double& length, int& error, spurt::vec3& start_dir, int nmax=10000)
{
    std::vector<spurt::vec3> dummy;
    return eigen_flow_map(field, seed, dx, length, error, start_dir, dummy, nmax);
}

template< typename Field >
spurt::vec3 eigen_flow_map(const Field& field, const spurt::vec3& seed, double dx, double& length, int& error,
                          std::vector<spurt::vec3>& steps, int nmax=10000)
{
    spurt::vec3 dummy(0, 0, 0);
    return eigen_flow_map(field, seed, dx, length, error, dummy, steps, nmax);
}

template< typename Field >
spurt::vec3 eigen_flow_map(const Field& field, const spurt::vec3& seed, double dx, double& length, int& error, int nmax=10000)
{
    spurt::vec3 dummy(0, 0, 0);
    return eigen_flow_map(field, seed, dx, length, error, dummy, nmax);
}

template< typename Field >
double ftle(int n, const Field& flowmap, double length)
{
    spurt::mat3 M(flowmap.derivative(flowmap.grid()(n)));

    spurt::mat3 T(M.transpose());
    M *= T;
    spurt::vec3 evals;
    spurt::mat3 evecs;
    sym_eigensystem(evals, evecs, M);
    return 1. / length*log(evals[0]);
}

template< typename Data, typename FlowMap >
spurt::vec2 eigenftle(int n, const Data& data, const FlowMap flowmaps[2], double length)
{
    auto coord = flowmaps[0].grid()(n);
    spurt::mat3 M[2];
    spurt::vec3 dmap[2];
    OrientationChecker<Data> checker(data, flowmaps[0].grid()(coord));
    for (int r = 0 ; r < 3 ; ++r) {
        central_diff(dmap, flowmaps, checker, coord, r);
        // std::cerr << "|dmap|(" << r << ") = " << spurt::norm(dmap[0]) << ", " << spurt::norm(dmap[1]) << std::endl;
        for (int c = 0 ; c < 3 ; ++c) {
            M[0](r, c) = dmap[0][c];
            M[1](r, c) = dmap[1][c];
        }
    }

    spurt::mat3 T[2] = { transpose(M[0]), transpose(M[1]) };
    double ftle[2] = { 0, 0 };
    spurt::vec3 evals;
    spurt::mat3 evecs;
    for (int i = 0 ; i < 2 ; ++i) {
        M[i] *= T[i];
        spurt::sym_eigensystem(evals, evecs, M[i]);
        ftle[i] = 1. / length * log(evals[0]);
    }
    spurt::vec2 ans(std::max(ftle[0], ftle[1]), std::min(ftle[0], ftle[1]));
    return ans;
}
}

#endif
