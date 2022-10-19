#include <poincare/map.hpp>

int xavier::map_debug::verbose_level = 0;
std::vector< nvis::vec2 > xavier::map_debug::jacobian_sample_pos;
std::vector< nvis::vec2 > xavier::map_debug::jacobian_sample_vals;
double xavier::map_debug::jacobian_timer;
