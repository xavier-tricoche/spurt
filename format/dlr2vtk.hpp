#ifndef __XAVIER_DLR_2_VTK_HPP__
#define __XAVIER_DLR_2_VTK_HPP__

// STL
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
// Boost
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
// VTK
#include <vtkCellArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetWriter.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
// nvis
#include <math/types.hpp>

#include "dlr_reader.hpp"
#include <vtk/vtk_utils.hpp>

namespace spurt {


class dlr2vtk_converter {

typedef fvec3                    point_type;
typedef bounding_box<point_type>  box_type;
typedef dlr_reader::cell_type   cell_type;
typedef std::pair<cell_type, long>      cell_entry;

typedef std::pair< std::vector<double>, std::string > pair_type;

// helper functions
struct cell_array_helper;

vtkPoints* vtk_points(const std::vector<point_type>&);

vtkUnstructuredGrid* vtk_cellarray(const cell_array_helper&);

public:

    dlr2vtk_converter() : m_boundary_only(false),
        m_mesh_only(false), m_verbose(false),
        m_comp_name("velocity") {}

    dlr2vtk_converter(const std::string& gridname="", const std::string& valname="")
        : m_grid_name(gridname), m_value_name(valname), m_attributes(),
          m_boundary_only(false),
          m_mesh_only(false), m_verbose(false), m_comp_name("velocity") {}

    void set_filenames(const std::string& gridname, const std::string& valname) {
        m_grid_name = gridname;
        m_value_name = valname;
    }

    void set_grid_name(const std::string& gridname) {
        m_grid_name = gridname;
    }

    void set_value_name(const std::string& valname) {
        m_value_name = valname;
    }

    void set_attributes(const std::vector<std::string>& attributes) {
        m_attributes = attributes;
    }

    void set_mesh_only(bool do_meshonly=true) {
        m_mesh_only = do_meshonly;
    }

    void set_boundary_only(bool do_boundary=true) {
        m_boundary_only = do_boundary;
    }

    void set_verbose(bool _verbose) {
        m_verbose = _verbose;
    }

    void import();

    VTK_SMART(vtkUnstructuredGrid) get_vtk_dataset() { return m_vtk_dataset; }

private:
    std::string m_grid_name, m_value_name, m_comp_name;
    std::vector<std::string> m_attributes;
    VTK_SMART(vtkUnstructuredGrid) m_vtk_dataset;
    bool m_boundary_only, m_mesh_only;
    bool m_verbose;
    std::vector< pair_type > m_vector_fields;
    std::vector< pair_type > m_scalar_fields;

    // helper functions
    struct cell_array_helper {
        cell_array_helper(const std::vector<long>& ids,
                          const std::vector<cell_entry>& types)
            : m_ids(ids), m_types(types) {}

        long start(const long& cell_id) const {
            return m_types[cell_id].second;
        }

        long end(const long& cell_id) const {
            return m_types[cell_id+1].second;
        }

        cell_type type(const long& cell_id) const {
            return m_types[cell_id].first;
        }

        long vertex(const long& loc_id, const long& cell_id) const {
            return m_ids[this->start(cell_id) + loc_id];
        }

        long vertex(const long& id) const {
            return m_ids[id];
        }

        long size(const long& cell_id) const {
            return this->end(cell_id) - this->start(cell_id);
        }

        long nb_cells() const {
            return m_types.size()-1;
        }

        const std::vector<long>&       m_ids;
        const std::vector<cell_entry>& m_types;
    };
};


void dlr2vtk_converter::import()
{
    // sanity check
    if (m_grid_name.empty()) {
        throw std::runtime_error("Missing grid file");
    }
    if (m_value_name.empty() && m_attributes.size())
        throw std::runtime_error("both value(s) and geometry were requested\n"
              "yet file containing values is missing");

    spurt::dlr_reader reader(m_grid_name, m_value_name);
    std::vector<point_type> points;
    std::vector<long> cell_indices;
    std::vector<cell_entry> cell_types;
    reader.read_mesh(m_boundary_only, points, cell_indices, cell_types, m_verbose);

    cell_array_helper helper(cell_indices, cell_types);

    if (!m_mesh_only) {
        for (int i=0; i<m_attributes.size() ; ++i) {
            std::string att = m_attributes[i];
            if (att[0] == '?' || att[0] == '*' || att[0] == '.') {
                if (att[1] == '_') att = att.substr(2);
                else att = att.substr(1);
                std::vector<double> val0, val1, val2;
                reader.read_data("x_" + att, val0, m_verbose);
                reader.read_data("y_" + att, val1, m_verbose);
                reader.read_data("z_" + att, val2, m_verbose);
                if (m_verbose)
                    std::cout << "val0.size()=" << val0.size() << '\n'
                        << "val1.size()=" << val1.size() << '\n'
                        << "val2.size()=" << val2.size() << '\n';
                assert(val0.size() == val1.size() && val1.size() == val2.size());
                m_vector_fields.push_back(pair_type());
                m_vector_fields.back().second = att;
                std::vector<double>& values = m_vector_fields.back().first;
                values.resize(3*val0.size());
                long k=0;
                for (long n=0; n<val0.size(); ++n) {
                    values[k++] = val0[n];
                    values[k++] = val1[n];
                    values[k++] = val2[n];
                }

                if (m_verbose) {
                    std::cout << "x_" << att << " range: "
                        << *std::min_element(val0.begin(), val0.end()) << " - "
                        << *std::max_element(val0.begin(), val0.end()) << '\n';

                    std::cout << "y_" << att << " range: "
                        << *std::min_element(val1.begin(), val1.end()) << " - "
                        << *std::max_element(val1.begin(), val1.end()) << '\n';

                    std::cout << "z_" << att << " range: "
                        << *std::min_element(val2.begin(), val2.end()) << " - "
                        << *std::max_element(val2.begin(), val2.end()) << '\n';
                }
            }
            else {
                m_scalar_fields.push_back(pair_type());
                std::vector<double>& values = m_scalar_fields.back().first;
                reader.read_data(att, values, m_verbose);
                m_scalar_fields.back().second = att;
            }
        }
    }

    box_type bounds;
    for (size_t i=0 ; i<points.size() ; ++i) {
        bounds.add(points[i]);
    }
    if (true || m_verbose) {
        std::cout << "there are " << points.size() << " points in input mesh\n"
                  << "associated bounding box: " << bounds << '\n';
    }

    if (m_verbose) std::cout << "about to create grid...\n";
    m_vtk_dataset = vtk_cellarray(helper);
    if (m_verbose) std::cout << "grid created. About to add points...\n";
    m_vtk_dataset->SetPoints(vtk_utils::make_vtkpoints(points));
    if (m_verbose) std::cout << "points added. About to add vectors...\n";
    for (int i=0; i<m_vector_fields.size(); ++i) {
        if (m_verbose)
            std::cout << "adding vector field #" << i << " of size " << m_vector_fields[i].first.size() << std::endl;
        vtk_utils::add_vectors_from_numbers<std::vector<double>, VTK_SMART(vtkUnstructuredGrid),3>
            (m_vtk_dataset, m_vector_fields[i].first, true, m_vector_fields[i].second);
    }
    if (m_verbose) std::cout << "vectors added. About to add scalars...\n";
    for (int i=0; i<m_scalar_fields.size(); ++i) {
        if (m_verbose) std::cout << "adding scalar field #" << i << "of size " << m_scalar_fields[i].first.size() << std::endl;
        vtk_utils::add_scalars(m_vtk_dataset, m_scalar_fields[i].first, true, m_scalar_fields[i].second);
    }
    if (m_verbose) std::cout << "scalars added\n";

    double* __bounds = m_vtk_dataset->GetBounds();
    box_type vtk_bounds(point_type(__bounds[0], __bounds[2], __bounds[4]),
                        point_type(__bounds[1], __bounds[3], __bounds[5]));

    if (m_verbose) {
        std::cout << "after conversion, exported grid has following properties:\n"
                  << "# points:     " << m_vtk_dataset->GetNumberOfPoints() << '\n'
                  << "# cells:      " << m_vtk_dataset->GetCells()->GetNumberOfCells() << '\n'
                  << "bounding box: " << vtk_bounds << '\n'
                  << "number of vectors: " << m_vtk_dataset->GetPointData()->GetScalars()->GetNumberOfTuples() << '\n'
                  << "number of scalars: " << m_vtk_dataset->GetPointData()->GetVectors()->GetNumberOfTuples() << '\n';
    }

}

vtkUnstructuredGrid*
dlr2vtk_converter::vtk_cellarray(const cell_array_helper& helper)
{
    if (m_verbose) std::cout << "in vtk_cellarray...\n";
    using namespace spurt;
    vtkCellArray* __cells = vtkCellArray::New();
    __cells->SetNumberOfCells(helper.nb_cells());
    for (long cell_id=0 ; cell_id<helper.nb_cells() ; ++cell_id) {
        __cells->InsertNextCell(helper.size(cell_id));
        for (long loc_id=0 ; loc_id<helper.size(cell_id) ; ++loc_id) {
            __cells->InsertCellPoint(helper.vertex(loc_id, cell_id));
        }
    }

    if (m_verbose) std::cout << "cell points inserted\n";

    vtkUnsignedCharArray* __types = vtkUnsignedCharArray::New();
    vtkIdTypeArray* __locations = vtkIdTypeArray::New();
    __types->SetNumberOfComponents(1);
    __types->SetNumberOfTuples(helper.nb_cells());
    __locations->SetNumberOfComponents(1);
    __locations->SetNumberOfTuples(helper.nb_cells());
    if (m_verbose) std::cout << "arrays initialized\n";
    for (long cell_id=0 ; cell_id<helper.nb_cells() ; ++cell_id) {
        unsigned char type_name;
        switch(helper.type(cell_id)) {
            case dlr_reader::TRIANGLE:
                type_name = VTK_TRIANGLE;
                break;
            case dlr_reader::QUADRILATERAL:
                type_name = VTK_QUAD;
                break;
            case dlr_reader::TETRAHEDRON:
                type_name = VTK_TETRA;
                break;
            case dlr_reader::HEXAHEDRON:
                type_name = VTK_HEXAHEDRON;
                break;
            case dlr_reader::PRISM:
                type_name = VTK_WEDGE;
                break;
            case dlr_reader::PYRAMID:
                type_name = VTK_PYRAMID;
                break;
            default:
                throw std::runtime_error("invalid cell type");
        }
        __types->SetValue(cell_id, type_name);
        __locations->SetValue(cell_id, helper.start(cell_id));
    }
    std::cout << "cell types defined\n";

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetCells(__types, __locations, __cells);
    std::cout << "grid cells set\n";
    __cells->Delete();
    std::cout << "cells pointer deleted\n";
    __types->Delete();
    std::cout << "types pointer deleted\n";
    __locations->Delete();
    std::cout << "locations pointer deleted\n";
    return grid;
}

} // namespace spurt


#endif
