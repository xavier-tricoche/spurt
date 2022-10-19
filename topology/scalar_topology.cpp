#include "topology/scalar_topology.hpp"
#include "image/nrrd_wrapper.hpp"

typedef xavier::topology::scalar::SmoothScalarField<double, 2> scalar_field_t;
typedef xavier::topology::scalar::planar_topology<double> topology_t;

int main(int argc, char* argv[]) {
    const std::string name = argv[1];
    
    Nrrd* nin = xavier::nrrd_utils::readNrrd(name);
    scalar_field_t field(nin);
    
    topology_t topology(field);
    std::vector<topology_t::critical_point> cps;
    topology.find_critical_points(cps, 50, 50);
    
    return 1;
}