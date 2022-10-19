#include "tokatopo.hpp"

int main(int argc, char* argv[])
{
    tokatopo::initialize(argc, argv);
    return tokatopo::compute_fixed_points();
}
