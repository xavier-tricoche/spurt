#ifndef __FTLE_HPP__
#define __FTLE_HPP__

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>
#include <math/matrix.hpp>
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
    OrientationChecker(const EigenvectorField& field, const nvis::vec3& x)
        : _field(field) {
        _ref_ev = _field.interpolate(x);
    }
    bool consistent(const nvis::vec3& x) const {
        nvis::vec3 ev = _field.interpolate(x);
        return (nvis::inner(_ref_ev, ev) > 0);
    }

    const EigenvectorField& _field;
    nvis::vec3 _ref_ev;
};

// helper function to compute flow map derivative on bi-directional flow
template < typename FlowMap, typename Checker >
void central_diff(nvis::vec3 ans[2], const FlowMap flowmaps[2],
                  const Checker& checker, const nvis::ivec3& coord, int dim)
{
    int i = coord[dim];
    int imax = flowmaps[0].grid().resolution()[dim];
    double h = flowmaps[0].grid().spacing()[dim];

    nvis::ivec3 ids[2] = { coord, coord }; // "left" and "right" neighbors
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
nvis::vec3 flow_map(const Field& field, const nvis::vec3& seed, double dx, double length)
{
    nvis::vec3 x = seed, lastx = seed;
    double l = 0, lastl = 0;
    while (l < length) {
        try {
            nvis::vec3 step = dx * field.interpolate(x);
            lastx = x;
            lastl = l;
            x += step;
            if (nvis::norm(step) < 1.0e-9) {
                return x;
            }
            l += nvis::norm(step);
        } catch (...) {
            return lastx;
        }
    }
    double u = (length - lastl) / (l - lastl);
    return (1. - u)*lastx + u*x;
}

template< typename Field >
nvis::vec3 eigen_flow_map(const Field& field, const nvis::vec3& seed, double dt, double& delta_T, int& error,
                          nvis::vec3& start_dir, std::vector<nvis::vec3>& steps, int nmax=10000)
{
    nvis::vec3 x = seed, lastx = seed;
    double t = 0, lastt = 0;

    steps.clear();
    steps.push_back(seed);

    try {
        // initialize integration direction
        nvis::vec3 ref_dir = field.interpolate(x);
        if (nvis::norm(start_dir) && nvis::inner(ref_dir, start_dir) < 0) {
            dt *= -1;
            ref_dir *= -1;
        }

        ref_dir /= nvis::norm(ref_dir); // ref_dir is always normalized

        int n=0;
        while (fabs(t) < delta_T && n<nmax) {
            nvis::vec3 step = dt * field.interpolate(x);
            // check orientation consistency
            double dot = nvis::inner(step, ref_dir);
            step *= spurt::sign(dot);
            start_dir = step;

            // std::cerr << "at " << x << ", t = " << t << ", step = " << step << std::endl;

            double step_length = nvis::norm(step);
            if (dot / step_length < 0.5) {
                // halve step size in case of sharp turn
                step *= 0.5;
                step_length *= 0.5;
            }
            if (step_length < 1.0e-9) {
                start_dir /= nvis::norm(start_dir);
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
        start_dir /= nvis::norm(start_dir);
        error = 0; // everything went fine
        double u = (delta_T - lastt) / (t - lastt);
        steps.push_back((1. - u)*lastx + u*x);
        return steps.back();
    } catch (std::runtime_error& e) {
        error = atoi(e.what());
        delta_T = lastt;
        start_dir /= nvis::norm(start_dir);
        return lastx;
    }
}

template< typename Field >
nvis::vec3 eigen_flow_map(const Field& field, const nvis::vec3& seed, double dx, double& length, int& error, nvis::vec3& start_dir, int nmax=10000)
{
    std::vector<nvis::vec3> dummy;
    return eigen_flow_map(field, seed, dx, length, error, start_dir, dummy, nmax);
}

template< typename Field >
nvis::vec3 eigen_flow_map(const Field& field, const nvis::vec3& seed, double dx, double& length, int& error,
                          std::vector<nvis::vec3>& steps, int nmax=10000)
{
    nvis::vec3 dummy(0, 0, 0);
    return eigen_flow_map(field, seed, dx, length, error, dummy, steps, nmax);
}

template< typename Field >
nvis::vec3 eigen_flow_map(const Field& field, const nvis::vec3& seed, double dx, double& length, int& error, int nmax=10000)
{
    nvis::vec3 dummy(0, 0, 0);
    return eigen_flow_map(field, seed, dx, length, error, dummy, nmax);
}

template< typename Field >
double ftle(int n, const Field& flowmap, double length)
{
    nvis::fixed_vector<nvis::vec3, 3> J(flowmap.derivative(flowmap.grid()(n)));
    spurt::mat3 M;
    for (int i = 0 ; i < 3 ; ++i) {
        for (int j = 0 ; j < 3 ; ++j) {
            M(i, j) = J[i][j];
        }
    }

    spurt::mat3 T(transpose(M));
    M *= T;
    std::vector<double> evals;
    std::vector<nvis::vec3> evecs;
    M.eigensystem(evals, evecs);
    if (evals.size()) {
        return 1. / length*log(*std::max_element(evals.begin(), evals.end()));
    } else {
        return 0.;
    }
}

template< typename Data, typename FlowMap >
nvis::vec2 eigenftle(int n, const Data& data, const FlowMap flowmaps[2], double length)
{
    nvis::ivec3 coord = flowmaps[0].grid()(n);
    spurt::mat3 M[2];
    nvis::vec3 dmap[2];
    OrientationChecker<Data> checker(data, flowmaps[0].grid()(coord));
    for (int r = 0 ; r < 3 ; ++r) {
        central_diff(dmap, flowmaps, checker, coord, r);
        // std::cerr << "|dmap|(" << r << ") = " << nvis::norm(dmap[0]) << ", " << nvis::norm(dmap[1]) << std::endl;
        for (int c = 0 ; c < 3 ; ++c) {
            M[0](r, c) = dmap[0][c];
            M[1](r, c) = dmap[1][c];
        }
    }

    spurt::mat3 T[2] = { transpose(M[0]), transpose(M[1]) };
    double ftle[2] = { 0, 0 };
    std::vector<double> evals;
    std::vector<nvis::vec3> evecs;
    for (int i = 0 ; i < 2 ; ++i) {
        M[i] *= T[i];
        M[i].eigensystem(evals, evecs);
        if (evals.size()) {
            ftle[i] = 1. / length * log(*std::max_element(evals.begin(), evals.end()));
        }
    }
    nvis::vec2 ans(std::max(ftle[0], ftle[1]), std::min(ftle[0], ftle[1]));
    return ans;
}
}

#endif
