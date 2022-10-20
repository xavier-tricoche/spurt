#ifndef __DATA_VTK_FIELD_HPP__
#define __DATA_VTK_FIELD_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <vtk/vtk_utils.hpp>
#include <vtkCellTreeLocator.h>
#include <vtkGenericCell.h>

namespace spurt {

class vtk_field {
    
enum AttributeType {
    NONE=-1,
    SCALAR,
    VECTOR,
    TENSOR,
    FIELD
};
    
template<typename T>
void interpolate(double pcoord[3], double* weights, T* res, vtkSmartPointer<vtkDataArray> values) const
{    
    vtkIdType id = m_cell->GetPointId(0);
    double tuple[9];
    values->GetTuple(id, tuple);
    
    for( int j=0; j<values->GetNumberOfComponents(); ++j )
        res[j] = tuple[j] * weights[0];

    for( int i=1; i<m_cell->GetNumberOfPoints(); ++i )
    {
        id = m_cell->GetPointId(i);
        values->GetTuple( id, tuple );
    
        for( int j=0; j<values->GetNumberOfComponents(); ++j )
            res[j] += tuple[j] * weights[i];
    }
}    
    
public:
    typedef double                         value_type;
    typedef nvis::fixed_vector<double, 3>  vector_type;
    typedef nvis::fixed_vector<double, 3>  point_type;
    typedef value_type                     scalar_type;
    typedef nvis::bounding_box<point_type> bounds_type;
    
    mutable size_t n_found, n_failed;
    
    vtk_field(const std::string& filename) {
        m_dataset = vtk_utils::readVTK(filename);
        m_locator = vtkSmartPointer<vtkCellTreeLocator>::New();
        m_locator->SetDataSet(m_dataset);
        m_locator->LazyEvaluationOff();
        m_locator->CacheCellBoundsOff();
        m_locator->BuildLocator();
        m_scalars = m_dataset->GetPointData()->GetScalars();
        m_vectors = m_dataset->GetPointData()->GetVectors();
        m_tensors = m_dataset->GetPointData()->GetTensors();
        m_cell = vtkSmartPointer<vtkGenericCell>::New();
        
        n_found = 0;
        n_failed = 0;
    }
    
    vtk_field(vtk_field& other) 
        : m_dataset(other.m_dataset), m_locator(other.m_locator), 
          m_scalars(other.m_scalars), m_vectors(other.m_vectors),
          m_tensors(other.m_tensors) {
        n_found = n_failed = 0;
        m_cell = vtkSmartPointer<vtkGenericCell>::New();
    }
    
    bounds_type bounds() const {
        double v[6];
        m_dataset->GetBounds(v);
        bounds_type b;
        b.min() = point_type(v[0], v[2], v[4]);
        b.max() = point_type(v[1], v[3], v[5]);
        return b;
    }
    
    template<typename Value_>
    bool operator()(const point_type& x, Value_& v, AttributeType att=SCALAR) const {
        double tol2, pcoords[3], weights[8];
        vtkIdType cellid = m_locator->FindCell(const_cast<double *>(&x[0]), tol2, m_cell, pcoords, weights);
        if (cellid == -1) {
            // ++n_failed;
            // std::cout << "cell search failed at " << p << "\n";
            return false;
        }
        m_dataset->GetCell(cellid, m_cell);
        // std::cout << "cell search succeeded at " << p << '\n';
        switch (att) {
        case SCALAR: interpolate(pcoords, weights, &v[0], m_scalars); break;
        case VECTOR: interpolate(pcoords, weights, &v[0], m_vectors); break;
        case TENSOR: interpolate(pcoords, weights, &v[0], m_tensors); break;
        default: throw std::runtime_error("Unrecognized attribute type" + std::to_string(att));
        }
        
        // ++n_found;
        return true;
    }

    vector_type operator()(const point_type& x) const {
        vector_type v;
        if (!this->operator()(x, v, VECTOR)) {
            throw std::runtime_error("invalid position: interpolation failed");
        }
        return v;
    }
    
private:
    vtkSmartPointer<vtkDataSet>             m_dataset;
    vtkSmartPointer<vtkCellTreeLocator>     m_locator;
    vtkSmartPointer<vtkDataArray>           m_scalars;
    vtkSmartPointer<vtkDataArray>           m_vectors;
    vtkSmartPointer<vtkDataArray>           m_tensors;
    mutable vtkSmartPointer<vtkGenericCell> m_cell;
};

} // spurt


#endif