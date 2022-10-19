#include "newton.hpp"
#include "fixpoints.hpp"

std::vector<nvis::vec2>                         xavier::newton_steps;
std::vector<std::pair<nvis::vec2, nvis::vec2> > xavier::search_steps;
bool                                            xavier::record_newton_steps;
bool                                            xavier::record_search_steps;
