/*
*   This is a hacked and augmented version of tendFiber.c,
*   from the Teem library as described below.
*
*    X. Tricoche and C. Garth, Jan 2007
*/

/*
    Teem: Tools to process and visualize scientific data and images
    Copyright (C) 2006, 2005  Gordon Kindlmann
    Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    (LGPL) as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.
    The terms of redistributing and/or modifying this software also
    include exceptions to the LGPL that facilitate static linking.
    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _INTEGRATE_HPP_
#define _INTEGRATE_HPP_

#include "teem/ten.h"
#include "math/fixed_vector.hpp"
#include <map>
#include <vector>
#include <assert.h>

#define __QUIET__

namespace xavier {

namespace ten_interface {

// Fiber: basic data structure to store tracks after integration

typedef nvis::vec3 point3;

struct HalfFiber : public std::vector< point3 > {
    HalfFiber() : total_length(0) {}

    double total_length;
    int why_stopped;
};

struct Fiber {
    Fiber();

    void reset();

    HalfFiber fwd, bwd;
    std::vector< point3 > samples_fwd, samples_bwd;
    point3 seed, initial_dir;
    bool valid;
};

std::ostream& operator<<(std::ostream& os, const Fiber& f)
{
    os << "FIBER: \n"
    << "\t valid = " << (f.valid ? "true" : "false") << '\n'
    << "\t initial_dir = " << f.initial_dir << '\n'
    << "\t forward = ";
    for (unsigned int i = 0 ; i < f.fwd.size() ; ++i) {
        os << i << ": " << f.samples_fwd[i] << " - ";
    }
    os << '\n'
    << "\t backward = ";
    for (unsigned int i = 0 ; i < f.bwd.size() ; ++i) {
        os << i << ": " << f.samples_bwd[i] << " - ";
    }
    os << '\n';

    return os;
}


// Chains are a sequential list of 3D points that support 1D spatial
// queries (probe). Piecewise linear interpolation is currently used to connect
// discrete points.
struct Chain : public std::vector< point3 > {
    Chain();
    Chain(unsigned int i0, unsigned int i1, const double* buffer);

    void probe(std::vector< point3 >& hf, const std::vector< double >& ts);

    std::vector< double > l;
};

// TrackingParam carries all the information relevant to the integration
// of tracks in a tensor volume. That includes the tensor volume itself,
// stop criteria, and integration control parameters.
struct TrackingParam {
    TrackingParam();
    ~TrackingParam();

    std::vector< point3 > seeds;
    double step;

    // Several stop criteria are provided corresponding to a compound
    // set of requirements on the integration. This follows the
    // approach used in tenFiber.c .
    // 0: aniso
    // 1: length
    // 2: nb steps
    std::vector< int > stop_criteria;

    // tenFiberIntg[Euler|Midpoint|RK4]
    int intg_method;

    // eigendirection to follow during integration
    // 0: major, 1: medium, 2:
    // NB: You should be using 0 for now
    int eigen;

    // tenAniso_[Cl1-2|Cp1-2|Ca1-2|Cs1-2|FA|Skew|Mode|Tr|eval0-2]
    int aniso_type;
    double min_aniso;
    double max_length;
    unsigned int max_nb_steps;
    double min_conf;
    double min_radius;

    Nrrd *dti;
};

struct context {
    context() : ptr(0) {}
    context(const TrackingParam& param);
    bool setup(const TrackingParam& param);
    ~context() {
        tenFiberContextNix(ptr);
        delete[] buffer;
    }

    tenFiberContext* ptr;
    unsigned int halfBuffLen;
    double *buffer;
};


void error_display(const std::string& who, const std::string& what,
                   const char* _type);

int integrate(const TrackingParam& param, std::vector< Fiber >& fibers,
              std::vector< double >& requested_lengths);

void integrate(Fiber& fiber, context& ctx);
void integrate_and_sample(Fiber& fiber, context& ctx,
                          const std::vector< double >& requested_lengths,
                          const bool save_curves = true);
};
};

#endif
