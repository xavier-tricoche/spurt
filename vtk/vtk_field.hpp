#ifndef __VTK_FIELD_HPP__
#define __VTK_FIELD_HPP__

#include "vtk_utils.hpp"

#include "vtkCell.h"
#include "vtkCellTreeLocator.h"
#include "vtkDataArray.h"
#include "vtkGenericCell.h"
#include "vtkIdList.h"
#include "vtkUnstructuredGrid.h"

#include <Eigen/Core>
#include <math/fixed_vector.hpp>

#include <string>


namespace spurt {

template<typename Scalar_>
class vtk_field {

template<typename T, typename NotUsed_> struct data_traits;

public:    
    typedef Scalar_ scalar_type;
    typedef fixed_vector<scalar_type, 3>     point_type;
    typedef fixed_vector<scalar_type, 3>     vector_type;
    typedef Eigen::Matrix<scalar_type, 3, 3> tensor_type;
    typedef nvis::bounding_box<point_type>   bounds_type;
    typedef vtk_field<scalar_type>           self_type;

    vtk_field(const std::string& filename, bool use_vectors=true, 
              bool use_scalars=false, bool use_tensors=false) 
        : m_has_scalars(use_scalars), m_has_vectors(use_vectors), m_has_tensors(use_tensors) {
        m_dataset = vtk_utils::readVTK(filename);
        
        if (m_has_scalars) {
            m_scalars = m_dataset->GetPointData()->GetScalars();
            if (!m_scalars) m_has_scalars = false;
        }
        if (m_has_vectors) {
            m_vectors = m_dataset->GetPointData()->GetVectors();
            if (!m_vectors) m_has_vectors = false;
        }
        if (m_has_tensors) {
            m_tensors = m_dataset->GetPointData()->GetTensors();
            if (!m_tensors) m_has_tensors = false;
        }
        
        // do we need a cell locator?
        if (vtkUnstructuredGrid::SafeDownCast(m_dataset)|| 
            vtkStructuredGrid::SafeDownCast(m_dataset)) {
            // the answer is yes
            m_locator = vtkCellTreeLocator::New();
            m_locator->SetDataSet(m_dataset);
            m_locator->SetCacheCellBounds(1);
            m_locator->AutomaticOn();
            m_locator->BuildLocator();
            m_has_locator = true;
        }
        else m_has_locator = false;
    }
    
    int locate(const point_type& x, VTK_SMART(vtkGenericCell) cell, double* weights) const {
        double pcoords[3];
        double tol2=1.0e-10;
        double xx[3] = {x[0], x[1], x[2]};
        return m_locator->FindCell(xx, tol2, cell, pcoords, weights);
    }
    
    template<typename T>
    T interpolate(const point_type& x) const {
        if (m_has_locator) {
            VTK_CREATE(vtkGenericCell, cell);
            double weights[8];
            int cell_id = locate(x, cell, weights);
            if (cell_id < 0) throw std::runtime_error("invalid position");
            T r;
            _M_interpolate(r, cell, weights);
            return r;
        }
        else {
            vtkCell* cell;
            int id, subid;
            double tol2 = 1.0e-10;
            double weights[8];
            double pcoords[3];
            double xx[3] = {x[0], x[1], x[2]};
            int cell_id = m_dataset->FindCell(xx, cell, id, tol2, subid, pcoords, weights);
            T r;
            _M_interpolate(r, cell_id, weights);
            return r;
        }
    }
    
    scalar_type interpolate_scalar(const point_type& x) const {
        if (!m_has_scalars) throw std::runtime_error("no scalar data available");
        return interpolate<scalar_type>(x);
    }

    vector_type interpolate_vector(const point_type& x) const {
        if (!m_has_vectors) throw std::runtime_error("no vector data available");
        return interpolate<vector_type>(x);
    }
    
    vector_type interpolate_tensor(const point_type& x) const {
        if (!m_has_tensors) throw std::runtime_error("no vector data available");
        return interpolate<tensor_type>(x);
    }
    
    void operator()(const point_type& x, vector_type& dxdt, scalar_type t=0) const {
        dxdt = interpolate_vector(x);
    }
    
    bounds_type bounds() const {
        double v[6];
        m_dataset->GetBounds(v);
        bounds_type b;
        b.min() = point_type(v[0], v[2], v[4]);
        b.max() = point_type(v[1], v[3], v[5]);
        return b;
    }
    
    
private:
    
    // convenient 
    template<typename T, typename NotUsed_ = void>
    struct data_traits {};
    
    template<typename Value_>
    void _M_interpolate(Value_& out, VTK_SMART(vtkGenericCell) cell,
                        double* weights) const {            
        VTK_SMART(vtkIdList) ids = cell->GetPointIds();
        _M_do_interpolate(out, ids, weights);
    }
    
    template<typename Value_>
    void _M_interpolate(Value_& out, int cellid, double* weights) const {            
        VTK_SMART(vtkIdList) ids;
         m_dataset->GetCellPoints(cellid, ids);
        _M_do_interpolate(out, ids, weights);
    }
    
    template<typename Value_>
    void _M_do_interpolate(Value_& out, VTK_SMART(vtkIdList) ids, double* weights) const {
        typedef Value_ value_type;
        typedef data_traits<Value_> value_traits;

        value_traits::set_to_zero(out);
        value_type val;
        double _val[9];
        for (int i=0; i<ids->GetNumberOfIds(); ++i) {
            value_traits::get_array(*this)->GetTuple(i, _val);
            out += weights[i]*val;
        }
    }
    
    template<typename NotUsed_> 
    struct data_traits<scalar_type, NotUsed_> {
        static void set_to_zero(scalar_type& s) { s=0; }
        static const VTK_SMART(vtkDataArray) get_array(const self_type& f) { return f.m_scalars; }
        static scalar_type* pointer(scalar_type& s) { return &s; }
    };

    template<typename NotUsed_> 
    struct data_traits<vector_type, NotUsed_> {
        static void set_to_zero(vector_type& v) { v[0]=v[1]=v[2]=0; }
        static const VTK_SMART(vtkDataArray) get_array(const self_type& f) { return f.m_vectors; }
        static scalar_type* pointer(vector_type& v) { return v.data(); }
    };

    template<typename NotUsed_> 
    struct data_traits<tensor_type, NotUsed_> {
        static void set_to_zero(tensor_type& t) { t=tensor_type::Zero(); }
        static const VTK_SMART(vtkDataArray) get_array(const self_type& f) { return f.m_tensors; }
        static scalar_type* pointer(tensor_type& t) { return &t(0,0); }
    };
    
    
    VTK_SMART(vtkDataSet) m_dataset;
    VTK_SMART(vtkCellTreeLocator) m_locator;
    VTK_SMART(vtkDataArray) m_scalars, m_vectors, m_tensors;
    bool m_has_locator, m_has_scalars, m_has_vectors, m_has_tensors;
};

} // namespace spurt


#endif
