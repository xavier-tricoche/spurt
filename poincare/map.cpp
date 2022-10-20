#include <poincare/map.hpp>

int spurt::map_debug::verbose_level = 0;
std::vector< nvis::vec2 > spurt::map_debug::jacobian_sample_pos;
std::vector< nvis::vec2 > spurt::map_debug::jacobian_sample_vals;
double spurt::map_debug::jacobian_timer;
