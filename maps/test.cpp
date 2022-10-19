#include "triangulation.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

typedef xavier::triangulation<double>       mesh_type;
typedef mesh_type::triangle_type            ids_type;

// export triangulation to file
void write(const mesh_type& tri,
           const std::string& name)
{
    std::fstream vtk(name.c_str(), std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << "triangulation test\n"
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";
        
    unsigned int nb_pts = tri.get_nb_vertices();
    vtk << "POINTS " << nb_pts << " float\n";
    // create a VTK point set with sampled points
    for (int i = 0 ; i < nb_pts ; ++i) {
        nvis::vec2 x = tri.get_vertex(i);
        vtk << x[0] << " " << x[1] << " 0.\n";
    }
    
    unsigned int ncells = tri.get_nb_triangles();
    vtk << "CELLS " << ncells << " " << 4*ncells << '\n';
    for (unsigned int n = 0 ; n < ncells ; ++n) {
        const ids_type ids = tri.get_triangle_vertices(n);
        vtk << "3 " << ids[0] << " " << ids[1] << " " << ids[2] << '\n';
    }
    
    vtk << "CELL_TYPES " << ncells << '\n';
    for (unsigned int i = 0 ; i < ncells ; ++i) {
        vtk << "5\n";
    }
    
    vtk << "POINT_DATA " << nb_pts << '\n'
        << "SCALARS radius float\n"
        << "LOOKUP_TABLE default\n";
    for (int i = 0 ; i < nb_pts ; ++i) {
        vtk << tri.get_data(i) << "\n";
    }
    vtk.close();
}

// triangulate a bunch of random points in unit square
int main(int argc, char* argv[])
{
    int nb_pts = atoi(argv[1]);
    
    std::vector<nvis::vec2> verts(nb_pts);
    std::vector<double> data(nb_pts);
    
    srand48(time(0));
    const nvis::vec2 center(0.5, 0.5);
    for (unsigned int i = 0 ; i < nb_pts ; ++i) {
        verts[i] = nvis::vec2(drand48(), drand48());
        data[i] = nvis::norm(verts[i] - center);
    }
    
    nvis::bbox2 bounds(nvis::vec2(0, 0), nvis::vec2(1, 1));
    mesh_type mesh(bounds, 0.);
    std::vector<mesh_type::index_type> modified;
    mesh.insert_points(verts, data, modified);
    
    write(mesh, argv[2]);
    
    return 0;
}





