#include "newton.hpp"
#include "fixpoints.hpp"

std::vector<nvis::vec2>                         spurt::newton_steps;
std::vector<std::pair<nvis::vec2, nvis::vec2> > spurt::search_steps;
bool                                            spurt::record_newton_steps;
bool                                            spurt::record_search_steps;
