#include <format/ncio.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <math/types.hpp>

int main(int argc, const char* argv[])
{
    std::string vtkname(argv[1]);
    size_t found = vtkname.rfind("nc");
    if (found != std::string::npos) {
        vtkname.replace(found, 2, "vtk");
    } else {
        std::cerr << "extension of " << argv[1] << " not recognized\n";
        return -1;
    }
    
    std::vector<double> vertices, velocity;
    std::vector<int> triangles;
    
    ncio::get(argv[1], vertices, "vertices");
    ncio::get(argv[1], triangles, "triangles");
    ncio::get(argv[1], velocity, "velocity");
    
    unsigned int npts = vertices.size() / 2;
    unsigned int ntri = triangles.size() / 3;
    unsigned int nvel = velocity.size() / 2;
    
    std::cerr << npts << " points, " << ntri << " triangles, " << nvel << " velocity vectors\n";
    
    std::fstream vtk(vtkname.c_str(), std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << "converted from NetCDF format file " << argv[1] << '\n'
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << npts << " float\n";
    for (int i = 0 ; i < npts ; ++i) {
        vtk << vertices[2*i] << " " << vertices[2*i+1] << " 0\n";
    }
    vtk << "CELLS " << ntri << " " << 4*ntri << '\n';
    for (int i = 0 ; i < ntri ; ++i) {
        vtk << "3";
        for (int j = 0 ; j < 3 ; ++j) {
            vtk << " " << triangles[3*i+j];
        }
        vtk << '\n';
    }
    vtk << "CELL_TYPES " << ntri << '\n';
    for (int i = 0 ; i < ntri ; ++i) {
        vtk << "5\n";
    }
    vtk <<  "POINT_DATA " << npts << '\n';
    vtk << "VECTORS velocity float\n";
    for (int i = 0 ; i < npts ; ++i) {
        vtk << velocity[2*i] << " " << velocity[2*i+1] << " 0\n";
    }
    vtk.close();
    
    return 0;
}





