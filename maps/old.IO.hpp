#ifndef __XAVIER_IO_HPP__
#define __XAVIER_IO_HPP__

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>


namespace spurt {

template<typename T>
inline void export_VTK(const T& mesh, const std::string& filename,
                       const std::string& comment, bool geometry_only = false, bool binary = false);
                       
template<typename T, typename Func>
inline void export_VTK(const T& mesh, const std::string& filename,
                       const std::string& comment, const Func& value_func, bool select = false);
                       
template<typename T, typename Func>
inline void export_submesh_VTK(const T& mesh, const std::string& filename,
                               const std::string& comment, const Func& value_func,
                               const std::vector<unsigned int>& ids, bool select = false);
}

template<typename T>
inline void spurt::export_VTK(const T& mesh, const std::string& filename,
                               const std::string& comment, bool geometry_only, bool binary)
{
    typedef T                                               triangulation_type;
    typedef typename triangulation_type::triangle_type      triangle_type;
    typedef typename triangulation_type::index_type         index_type;
    typedef typename triangulation_type::point_type         point_type;
    typedef typename triangulation_type::data_type          data_type;
    
    std::fstream vtk(filename.c_str(), std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << comment << '\n'
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";
        
    size_t nb_pts = mesh.get_nb_vertices();
    vtk << "POINTS " << nb_pts << " float\n";
    for (size_t i = 0 ; i < nb_pts ; ++i) {
        const point_type& p = mesh.get_vertex(i);
        vtk << p[0] << " " << p[1] << " 0.\n";
    }
    
    size_t ncells = mesh.get_nb_triangles();
    vtk << "CELLS " << ncells << " " << 4*ncells << '\n';
    for (size_t i = 0 ; i < ncells ; ++i) {
        const triangle_type& tri = mesh.get_triangle_vertices(i);
        vtk << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }
    
    vtk << "CELL_TYPES " << ncells << '\n';
    for (size_t i = 0 ; i < ncells ; ++i) {
        vtk << "5\n";
    }
    
    if (!geometry_only) {
        vtk << "POINT_DATA " << nb_pts << '\n'
            << "SCALARS q float\n"
            << "LOOKUP_TABLE default\n";
        for (size_t i = 0 ; i < nb_pts ; ++i) {
            if (geometry_only) {
                vtk << "0\n";
            } else {
                const data_type& v = mesh.get_data(i);
                vtk << v.scalar_value() << '\n';
            }
        }
    }
    
    vtk.close();
}

template<typename T, typename Func>
inline void spurt::export_VTK(const T& mesh, const std::string& filename,
                               const std::string& comment, const Func& value_func, bool select)
{
    typedef T                                               triangulation_type;
    typedef typename triangulation_type::triangle_type      triangle_type;
    typedef typename triangulation_type::index_type         index_type;
    typedef typename triangulation_type::point_type         point_type;
    typedef typename triangulation_type::data_type          data_type;
    
    
    size_t nb_pts = mesh.get_nb_vertices();
    size_t ncells = mesh.get_nb_triangles();
    std::vector<int> old2new(nb_pts);
    
    std::vector<bool> include_vertex(nb_pts, true);
    std::vector<bool> include_cell(ncells, true);
    int counter = 0;
    if (select) {
        for (int i = 0 ; i < include_vertex.size() ; ++i) {
            include_vertex[i] = value_func.is_valid(mesh.get_data(i));
            if (include_vertex[i]) {
                old2new[i] = counter++;
            } else {
                --nb_pts;
            }
        }
        
        for (int i = 0 ; i < include_cell.size() ; ++i) {
            const triangle_type& tri = mesh.get_triangle_vertices(i);
            include_cell[i] = (include_vertex[tri[0]] &&
                               include_vertex[tri[1]] &&
                               include_vertex[tri[2]]);
            if (!include_cell[i]) {
                --ncells;
            }
        }
    }
    
    std::fstream vtk(filename.c_str(), std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << comment << '\n'
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";
        
    vtk << "POINTS " << nb_pts << " float\n";
    for (size_t i = 0 ; i < include_vertex.size() ; ++i) {
        if (!include_vertex[i]) {
            continue;
        }
        const point_type& p = mesh.get_vertex(i);
        vtk << p[0] << " " << p[1] << " 0.\n";
    }
    
    vtk << "CELLS " << ncells << " " << 4*ncells << '\n';
    for (size_t i = 0 ; i < include_cell.size() ; ++i) {
        if (!include_cell[i]) {
            continue;
        }
        const triangle_type& tri = mesh.get_triangle_vertices(i);
        vtk << "3 " << old2new[tri[0]] << " " << old2new[tri[1]] << " " << old2new[tri[2]] << "\n";
    }
    
    vtk << "CELL_TYPES " << ncells << '\n';
    for (size_t i = 0 ; i < ncells ; ++i) {
        vtk << "5\n";
    }
    
    int order = value_func.order();
    vtk << "POINT_DATA " << nb_pts << '\n';
    if (!order) {
        vtk << "SCALARS " << value_func.name() << " float\n"
            << "LOOKUP_TABLE default\n";
    } else if (order == 1) {
        vtk << "VECTORS " << value_func.name() << " float\n";
    } else if (order == 2) {
        vtk << "TENSORS " << value_func.name() << " float\n";
    } else {
        assert(false);
    }
    for (size_t i = 0 ; i < include_vertex.size() ; ++i) {
        if (!include_vertex[i]) {
            continue;
        }
        const data_type& v = mesh.get_data(i);
        vtk << value_func.value_string(v) << '\n';
    }
    
    vtk.close();
}

template<typename T, typename Func>
inline void spurt::export_submesh_VTK(const T& mesh, const std::string& filename,
                                       const std::string& comment, const Func& value_func,
                                       const std::vector<unsigned int>& ids, bool select)
{
    typedef T                                               triangulation_type;
    typedef typename triangulation_type::triangle_type      triangle_type;
    typedef typename triangulation_type::index_type         index_type;
    typedef typename triangulation_type::point_type         point_type;
    typedef typename triangulation_type::data_type          data_type;
    
    // old to new and back
    std::map<unsigned int, unsigned int> old2new;
    std::map<unsigned int, unsigned int>::iterator it;
    std::vector<unsigned int> new_to_old;
    
    // loop over all valid triangles to identify included vertices
    for (int i = 0 ; i < ids.size() ; ++i) {
        triangle_type tri = mesh.get_triangle_vertices(ids[i]);
        for (int j = 0 ; j < 3 ; ++j) {
            it = old2new.find(tri[j]);
            if (it == old2new.end()) {
                nvis::vec2 x = mesh.get_vertex(tri[j]);
                new_to_old.push_back(tri[j]);
                old2new[tri[j]] = new_to_old.size() - 1;
            }
        }
    }
    size_t nb_pts = new_to_old.size();
    size_t ncells = ids.size();
    
    std::fstream vtk(filename.c_str(), std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << comment << '\n'
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";
        
    vtk << "POINTS " << nb_pts << " float\n";
    for (size_t i = 0 ; i < new_to_old.size() ; ++i) {
        const point_type& p = mesh.get_vertex(new_to_old[i]);
        vtk << p[0] << " " << p[1] << " 0.\n";
    }
    
    vtk << "CELLS " << ncells << " " << 4*ncells << '\n';
    for (size_t i = 0 ; i < ids.size() ; ++i) {
        const triangle_type& tri = mesh.get_triangle_vertices(ids[i]);
        vtk << "3 " << old2new[tri[0]] << " " << old2new[tri[1]] << " " << old2new[tri[2]] << "\n";
    }
    
    vtk << "CELL_TYPES " << ncells << '\n';
    for (size_t i = 0 ; i < ncells ; ++i) {
        vtk << "5\n";
    }
    
    int order = value_func.order();
    vtk << "POINT_DATA " << nb_pts << '\n';
    if (!order) {
        vtk << "SCALARS " << value_func.name() << " float\n"
            << "LOOKUP_TABLE default\n";
    } else if (order == 1) {
        vtk << "VECTORS " << value_func.name() << " float\n";
    } else if (order == 2) {
        vtk << "TENSORS " << value_func.name() << " float\n";
    } else {
        assert(false);
    }
    for (size_t i = 0 ; i < new_to_old.size() ; ++i) {
        const data_type& v = mesh.get_data(new_to_old[i]);
        vtk << value_func.value_string(v) << '\n';
    }
    
    vtk.close();
}






#endif


















