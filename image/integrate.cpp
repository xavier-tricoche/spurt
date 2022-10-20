/*
*   This is a hacked and augmented version of tendFiber.c,
*   from the Teem library as described below.
*
*    X. Tricoche and C. Garth, Jan - Oct 2007
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

#include "integrate.hpp"
#include <iostream>
#include "misc/progress.hpp"

using namespace spurt;
using namespace ten_interface;

typedef spurt::ProgressDisplay progress_t;

spurt::ten_interface::Fiber::Fiber()
        : fwd(), bwd(), valid(true)
{}

void spurt::ten_interface::Fiber::reset()
{
    fwd.clear();
    bwd.clear();
    fwd.total_length = 0;
    bwd.total_length = 0;

    initial_dir = nvis::vec3(0, 0, 0);

    valid = true;

    fwd.why_stopped = tenFiberStopUnknown;
    bwd.why_stopped = tenFiberStopUnknown;

    samples_fwd.clear();
    samples_bwd.clear();
}

spurt::ten_interface::Chain::Chain() {}

spurt::ten_interface::Chain::Chain(unsigned int i0, unsigned int i1, const double* buffer)
{
    this->resize((i1 > i0 ? i1 - i0 + 1 : i0 - i1 + 1));
    if (i1 > i0)
        for (unsigned int i = 0 ; i <= i1 - i0 ; i++) {
            (*this)[i][0] = buffer[3*(i+i0)];
            (*this)[i][1] = buffer[3*(i+i0)+1];
            (*this)[i][2] = buffer[3*(i+i0)+2];
        }
    else
        for (unsigned int i = 0 ; i <= i0 - i1 ; i++) {
            (*this)[i][0] = buffer[3*(i0-i)];
            (*this)[i][1] = buffer[3*(i0-i)+1];
            (*this)[i][2] = buffer[3*(i0-i)+2];
        }

    l.resize(size());
    l[0] = 0;
    for (unsigned int i = 1 ; i < size() ; i++) {
        l[i] = l[i-1] + nvis::norm((*this)[i] - (*this)[i-1]);
    }
}

void spurt::ten_interface::Chain::probe(std::vector< point3 >& hf, const std::vector< double >& ts)
{
    assert(ts.size() && ts[0] > 0);

    hf.clear();
    unsigned int ti = 0;
    unsigned int li = 1; // l[li-1]<ts[ti];
    while (ti < ts.size() && ts[ti] < l.back()) {
        for (; li < l.size() && l[li] < ts[ti] ; li++) {}
        if (li >= l.size()) break;

        double u = (ts[ti] - l[li-1]) / (l[li] - l[li-1]);
        hf.push_back((1 - u)*(*this)[li-1] + u*(*this)[li]);

        ti++;
    }
}

spurt::ten_interface::TrackingParam::TrackingParam()
        : seeds(), step(0.01), stop_criteria(), intg_method(tenFiberIntgRK4),
        eigen(tenFiberTypeEvec0),
        aniso_type(tenAniso_FA), min_aniso(0.5), max_length(10.0),
        max_nb_steps(500), min_conf(0.5), min_radius(100.), dti(0)
{
    stop_criteria.resize(5);
    stop_criteria[0] = 0;
    stop_criteria[1] = 1;
    stop_criteria[2] = 2;
    stop_criteria[3] = 3;
    stop_criteria[4] = 4;
}

spurt::ten_interface::TrackingParam::~TrackingParam()
{
    dti = 0;
}

void spurt::ten_interface::
error_display(const std::string& who, const std::string& what,
              const char* _type)
{
    char *err = biffGetDone(_type);
    std::cout << who << " reported: " << what << std::endl;
    std::cout << "internal error message: "
    << err << std::endl;
}

int spurt::ten_interface::
integrate(const TrackingParam& param, std::vector< Fiber >& fibers,
          std::vector< double >& requested_lengths)
{
    tenFiberContext *tfx;
    double start[3];
    int E;

    std::cout << "used integration lengths: " << std::endl;
    for (unsigned int i = 0 ; i < requested_lengths.size() ; i++)
        std::cout << requested_lengths[i] << "mm ";
    std::cout << std::endl;

//    double threshold = 0;

    tfx = tenFiberContextNew(param.dti);
    if (tfx->useIndexSpace)
        std::cout << "WE ARE USING INDEX SPACE" << std::endl;
    else
        std::cout << "WE ARE USING WORLD COORDINATES" << std::endl;

    if (!tfx)
        ten_interface::error_display("integrate", "failed to create fiber context", TEN);

    // assign stop criteria
    E = 0;
    for (unsigned int si = 0 ; si < param.stop_criteria.size() ; si++) {
        switch (param.stop_criteria[si]) {
        case 0:
            if (!E) E |= tenFiberStopSet(tfx, tenFiberStopAniso,
                                             param.aniso_type,
                                             param.min_aniso);
            break;
        case 1:
            if (!E) E |= tenFiberStopSet(tfx, tenFiberStopLength,
                                             param.max_length);
            break;
        case 2:
            if (!E) E |= tenFiberStopSet(tfx, tenFiberStopNumSteps,
                                             (int)param.max_nb_steps);
            break;
        case 3:
            if (!E) E |= tenFiberStopSet(tfx, tenFiberStopConfidence,
                                             param.min_conf);
            break;
        case 4:
            if (!E) E |= tenFiberStopSet(tfx, tenFiberStopRadius,
                                             param.min_radius);
            break;
        case tenFiberStopBounds:
        default:
            /* nothing to do */
            break;
        }
    }

    if (!E) E |= tenFiberTypeSet(tfx, param.eigen);  // considered eigenvector

    // hard-coded kernel for now: smoothing BC kernel
    double kparms[NRRD_KERNEL_PARMS_NUM];
    kparms[0] = 1.0;
    kparms[1] = 1.0;
    kparms[2] = 0.0;
    if (!E) E |= tenFiberKernelSet(tfx, nrrdKernelBCCubic, kparms);

    if (!E) E |= tenFiberIntgSet(tfx, param.intg_method);
    if (!E) E |= tenFiberParmSet(tfx, tenFiberParmStepSize, param.step); // step size
    if (!E) E |= tenFiberParmSet(tfx, tenFiberParmUseIndexSpace, AIR_FALSE);
    if (!E) E |= tenFiberUpdate(tfx);

    if (E)
        ten_interface::error_display("integrate", "trouble setting up ten fiber context", TEN);

    progress_t t0;
    t0.start();
    unsigned int nbfib = 0;
    unsigned int halfBuffLen = param.max_nb_steps;
    double *buff = new double[3*(2*halfBuffLen+1)];
    unsigned int id0, id1;

    unsigned int nfibers = param.seeds.size();
    fibers.resize(nfibers);

    unsigned int valid_dir = 0;

    std::cout << "number of fibers: " << nfibers << std::endl;

    for (unsigned int i = 0 ; i < nfibers ; i++) {
        start[0] = param.seeds[i][0];
        start[1] = param.seeds[i][1];
        start[2] = param.seeds[i][2];

        fibers[i].reset();

        id0 = id1 = halfBuffLen;

        // each half fiber should at least know about its seed point and
        // initial direction regardless the outcome of the numerical
        // integration
        fibers[i].seed = param.seeds[i];
        fibers[i].initial_dir = point3(0, 0, 0);

        if (tenFiberTraceSet(tfx, NULL, buff, halfBuffLen, &id0, &id1, start)) {
            ten_interface::error_display("integrate", "trouble with tracking", TEN);
            std::cout << "we got screwed" << std::endl;
            fibers[i].valid = false;
        }
        else {
            if (id0 < halfBuffLen && id0 < id1 && id1 > halfBuffLen) {
                ++nbfib;
                Chain c1(halfBuffLen, id1 - 1, buff);
                Chain c2(halfBuffLen, id0 + 1, buff);
                fibers[i].fwd.total_length = c1.l.back();
                fibers[i].bwd.total_length = c2.l.back();
                std::copy(c1.begin(), c1.end(), std::back_inserter(fibers[i].fwd));
                std::copy(c2.begin(), c2.end(), std::back_inserter(fibers[i].bwd));
                point3 dir = c1[1] - c2[1];
                dir *= 1 / norm(dir);
                fibers[i].initial_dir = dir;

                ++valid_dir;
                c1.probe(fibers[i].samples_fwd, requested_lengths);
                c2.probe(fibers[i].samples_bwd, requested_lengths);
                // check if we were able to reach at least the shortest requested length
                // in both directions
                if (!fibers[i].samples_fwd.size() || !fibers[i].samples_bwd.size()) {
                    fibers[i].valid = false;
                }
                // std::cout << "integration succeeded!" << std::endl;
            }
            else {
                fibers[i].valid = false;
            }
        }
    }

    std::cout << "we were able to assign " << valid_dir << " valid orientations "
    << "for a total of " << fibers.size() << " fibers" << std::endl;
    t0.end();
    std::cout << (float)nbfib / t0.wall_time() << " fibers per second" << std::endl;

    return 0;
}


spurt::ten_interface::context::
context(const TrackingParam& param)
{
    int E;
    halfBuffLen = param.max_nb_steps;
    buffer = new double[3*(2*halfBuffLen+1)];

    ptr = tenFiberContextNew(param.dti);

    if (!ptr) {
        ten_interface::error_display("integrate", "failed to create fiber context", TEN);
        exit(-1);
    }

    // assign stop criteria
    E = 0;
    for (unsigned int si = 0 ; si < param.stop_criteria.size() ; si++) {
        switch (param.stop_criteria[si]) {
        case 0:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopAniso,
                                             param.aniso_type,
                                             param.min_aniso);
            break;
        case 1:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopLength,
                                             param.max_length);
            break;
        case 2:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopNumSteps,
                                             (int)param.max_nb_steps);
            break;
        case 3:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopConfidence,
                                             param.min_conf);
            break;
        case 4:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopRadius,
                                             param.min_radius);
            break;
        case tenFiberStopBounds:
        default:
            /* nothing to do */
            break;
        }
    }

    if (!E) E |= tenFiberTypeSet(ptr, param.eigen);  // considered eigenvector

    // hard-coded kernel for now: smoothing BC kernel
    double kparms[NRRD_KERNEL_PARMS_NUM];
    kparms[0] = 1.0;
    kparms[1] = 1.0;
    kparms[2] = 0.0;
    if (!E) E |= tenFiberKernelSet(ptr, nrrdKernelBCCubic, kparms);
    if (!E) E |= tenFiberIntgSet(ptr, param.intg_method);
    if (!E) E |= tenFiberParmSet(ptr, tenFiberParmStepSize, param.step); // step size
    if (!E) E |= tenFiberParmSet(ptr, tenFiberParmUseIndexSpace, AIR_FALSE);
    if (!E) E |= tenFiberUpdate(ptr);

    if (E) {
        ten_interface::error_display("context constructor", "trouble setting up ten fiber context", TEN);
        exit(-1);
    }
}

bool spurt::ten_interface::context::
setup(const TrackingParam& param)
{
    int E;
    halfBuffLen = param.max_nb_steps;
    if (buffer) delete[] buffer;
    buffer = new double[3*(2*halfBuffLen+1)];

    if (ptr) tenFiberContextNix(ptr);
    ptr = tenFiberContextNew(param.dti);
    if (ptr->useIndexSpace)
        std::cout << "WE ARE USING INDEX SPACE RIGHT NOW" << std::endl;

    if (!ptr) {
        ten_interface::error_display("integrate", "failed to create fiber context", TEN);
        return false;
    }

    // assign stop criteria
    E = 0;
    for (unsigned int si = 0 ; si < param.stop_criteria.size() ; si++) {
        switch (param.stop_criteria[si]) {
        case 0:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopAniso,
                                             param.aniso_type,
                                             param.min_aniso);
            break;
        case 1:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopLength,
                                             param.max_length);
            break;
        case 2:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopNumSteps,
                                             (int)param.max_nb_steps);
            break;
        case 3:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopConfidence,
                                             param.min_conf);
            break;
        case 4:
            if (!E) E |= tenFiberStopSet(ptr, tenFiberStopRadius,
                                             param.min_radius);
            break;
        case tenFiberStopBounds:
        default:
            /* nothing to do */
            break;
        }
    }

    if (!E) E |= tenFiberTypeSet(ptr, param.eigen);  // considered eigenvector

    // hard-coded kernel for now: smoothing BC kernel
    double kparms[NRRD_KERNEL_PARMS_NUM];
    kparms[0] = 1.0;
    kparms[1] = 1.0;
    kparms[2] = 0.0;
    if (!E) E |= tenFiberKernelSet(ptr, nrrdKernelBCCubic, kparms);

    if (!E) E |= tenFiberIntgSet(ptr, param.intg_method);
    if (!E) E |= tenFiberParmSet(ptr, tenFiberParmStepSize, param.step); // step size
    if (!E) E |= tenFiberParmSet(ptr, tenFiberParmUseIndexSpace, AIR_FALSE);
    if (!E) E |= tenFiberUpdate(ptr);

    if (E) {
        ten_interface::error_display("integrate", "trouble setting up ten fiber context", TEN);
        return false;
    }

    return true;
}

void spurt::ten_interface::integrate(Fiber& fiber, context& ctx)
{
    double start[3] = { fiber.seed[0], fiber.seed[1], fiber.seed[2] };
    unsigned int id0, id1;

    fiber.reset();

    id0 = id1 = ctx.halfBuffLen;

    // each half fiber should at least know about its seed point and
    // initial direction regardless the outcome of the numerical
    // integration
    fiber.initial_dir = point3(0, 0, 0);

    if (tenFiberTraceSet(ctx.ptr, NULL, ctx.buffer, ctx.halfBuffLen, &id0, &id1, start)) {
        ten_interface::error_display("integrate", "trouble with tracking", TEN);
        std::cout << "we got screwed" << std::endl;
        fiber.valid = false;
    }
    else {
        if (id0 < ctx.halfBuffLen && id0 < id1 && id1 > ctx.halfBuffLen) {
            Chain c1(ctx.halfBuffLen, id1 - 1, ctx.buffer);
            Chain c2(ctx.halfBuffLen, id0 + 1, ctx.buffer);
            point3 dir = c1[1] - c2[1];
            dir *= 1 / norm(dir);
            fiber.initial_dir = dir;
            std::copy(c1.begin(), c1.end(), std::back_inserter(fiber.fwd));
            std::copy(c2.begin(), c2.end(), std::back_inserter(fiber.bwd));
            fiber.fwd.total_length = c1.l.back();
            fiber.bwd.total_length = c2.l.back();
            fiber.valid = true;
        }
        else {
            fiber.valid = false;
        }
        fiber.fwd.why_stopped = ctx.ptr->whyStop[1];
        fiber.bwd.why_stopped = ctx.ptr->whyStop[0];
    }
}

void spurt::ten_interface::integrate_and_sample(Fiber& fiber, context& ctx,
        const std::vector< double >& requested_lengths, const bool save_curves)
{
    double start[3] = { fiber.seed[0], fiber.seed[1], fiber.seed[2] };
    unsigned int id0, id1;

    // fiber.reset();

    id0 = id1 = ctx.halfBuffLen;

    // each half fiber should at least know about its seed point and
    // initial direction regardless the outcome of the numerical
    // integration
    fiber.initial_dir = point3(0, 0, 0);

    if (tenFiberTraceSet(ctx.ptr, NULL, ctx.buffer, ctx.halfBuffLen, &id0, &id1, start)) {
        ten_interface::error_display("integrate", "trouble with tracking", TEN);
        std::cout << "we got screwed" << std::endl;
    }
    else {
        if (id0 < ctx.halfBuffLen && id0 < id1 && id1 > ctx.halfBuffLen) {
            fiber.valid = true;

            Chain c1(ctx.halfBuffLen, id1 - 1, ctx.buffer);
            Chain c2(ctx.halfBuffLen, id0 + 1, ctx.buffer);
            // copy integral curves
            if (save_curves) {
                std::copy(c1.begin(), c1.end(), std::back_inserter(fiber.fwd));
                std::copy(c2.begin(), c2.end(), std::back_inserter(fiber.bwd));
            }
            // initial direction for orientation consistency
            point3 dir = c1[1] - c2[1];
            dir *= 1 / norm(dir);
            fiber.initial_dir = dir;
            // integration diagnostic
            fiber.fwd.why_stopped = ctx.ptr->whyStop[1];
            fiber.bwd.why_stopped = ctx.ptr->whyStop[0];
            // requested positions along curves
            c1.probe(fiber.samples_fwd, requested_lengths);
            c2.probe(fiber.samples_bwd, requested_lengths);
            // check if we were able to reach at least the shortest requested length
            // in both directions
            if (!fiber.samples_fwd.size() || !fiber.samples_bwd.size()) {
                fiber.valid = false;
            }
        }
        else {
            fiber.valid = false;
            fiber.fwd.why_stopped = ctx.ptr->whyStop[1];
            fiber.bwd.why_stopped = ctx.ptr->whyStop[0];
        }
    }
}
