#pragma once

#include <math/types.hpp>
#include <image/probe.hpp>
#include <image/eigen.hpp>
#include <vector>
#include <list>

namespace spurt {
namespace crease {
// this parameter controls the subdivision going on in the
// detection of crease points on provided edges.
// 0 value will lead to no subdivision
unsigned int nsub;

bool find_intersection(double* x0, double* x1, double* inter,
                       gage_interface::scalar_wrapper& gH_wrapper,
                       bool ridge);
bool find_intersection(const vec2& x0, const vec2& x1,
                       vec2& inter,
                       gage_interface::scalar_wrapper& gH_wrapper,
                       bool ridge);
                       
bool find_intersection(const vec2& x0, const vec2& x1,
                       vec2& inter,
                       gage_interface::scalar_wrapper& gH_wrapper,
                       bool ridge, double& strength);
                       
bool quick_search(double* x0, double* x1, double* inter,
                  gage_interface::scalar_wrapper& gH,
                  bool ridge);
bool quick_search(const vec2& x0, const vec2& x1,
                  vec2& inter,
                  gage_interface::scalar_wrapper& gH,
                  bool ridge);
                  
bool eigen(double* emin, double& lmin, double* H,
           bool ridge);
bool eigen(vec2& emin, double& lmin, const vec3& H,
           bool ridge);
};
};