#ifndef __VTK_INTERPOLATOR_HPP__
#define __VTK_INTERPOLATOR_HPP__

#include <VTK/vtk_utils.hpp>

#include <vtkGenericCell.h>
#include <vtkCellLocator.h>
#include <vtkCellTreeLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkImageData.h>
#include <vtkStructuredGrid.h>
#include <vtkPixel.h>

#include <cstdlib>

namespace {
    struct debug_info {
        bool verbose;
        std::ostringstream os;

        debug_info() : verbose(false), os() {}
        debug_info(bool v, const std::ostringstream& _os=std::ostringstream()) : verbose(v), os() {
            os << _os.str();
        }
    };

    template<typename S>
    struct is_double_valued : public std::false_type {};

    template<>
    struct is_double_valued<double> : public std::true_type {};

    template<typename DataSet>
    struct is_structured : public std::false_type {};
    template<>
    struct is_structured<vtkImageData> : public std::true_type {};
    template<>
    struct is_structured<vtkUniformGrid> : public std::true_type {};
    template<>
    struct is_structured<vtkRectilinearGrid> : public std::true_type {};
    template<>
    struct is_structured<vtkStructuredGrid> : public std::true_type {};
    template<>
    struct is_structured<boundaryAwareRectGrid> : public std::true_type {}; //BARG edit

    template<typename DataSet>
    struct needs_cell_locator : public std::false_type {};
    template<>
    struct needs_cell_locator<vtkUnstructuredGrid> : public std::true_type {};
    template<>
    struct needs_cell_locator<vtkPolyData> : public std::true_type {};
    template<>
    struct needs_cell_locator<vtkStructuredGrid> : public std::true_type {};


    template<typename Dataset, typename Enable=void >
    struct index_converter {};

    template<typename Dataset>
    struct index_converter< Dataset, typename std::enable_if<is_structured<Dataset>::value>::type > {
        index_converter(Dataset* data) {
            data->GetDimensions(m_dims);
            int dx = 1;
            int dy = m_dims[0];
            int dz = m_dims[0]*m_dims[1];
            m_deltas[0] = 0;
            m_deltas[1] = dx;
            m_deltas[2] = dy;
            m_deltas[3] = dy+dx;
            m_deltas[4] = dz;
            m_deltas[5] = dz+dx;
            m_deltas[6] = dz+dy;
            m_deltas[7] = dz+dy+dx;
        }

        int cellid_2_pointid(int cellid) const {
            std::div_t qr0 = std::div(cellid, m_dims[0]-1);
            std::div_t qr1 = std::div(qr0.quot, m_dims[1]-1);
            const int& i = qr0.rem;
            const int& j = qr1.rem;
            const int& k = qr1.quot;
            return i + m_dims[0]*(j + m_dims[1]*k);
        }

        void cellid_2_vertids(int cellid, std::vector<int>& verts) const {
            std::div_t qr0 = std::div(cellid, m_dims[0]-1);
            std::div_t qr1 = std::div(qr0.quot, m_dims[1]-1);
            const int& i = qr0.rem;
            const int& j = qr1.rem;
            const int& k = qr1.quot;
            verts.resize(8);
            verts[0] = i + m_dims[0]*(j + m_dims[1]*k);
            for (int i=1; i<8; ++i) {
                verts[i] = verts[0] + m_deltas[i];
            }
        }

        int delta(int i) const {
            return m_deltas[i];
        }

        int dim(int i) const {
            return m_dims[i];
        }

        int m_dims[3];
        int m_deltas[8];
    };

    template<typename V, typename S=typename V::value_type, typename Enable=void>
    struct find_cell_without_locator {};

    template<typename V, typename S>
    struct find_cell_without_locator<V, S,
        typename std::enable_if<is_double_valued<S>::value>::type> {
        template<typename DataSet>
        int find(const V& x, double* pcoord, double* weights, const DataSet* ds) {
            int subid;
            return const_cast<DataSet*>(ds)->FindCell(const_cast<S*>(&x[0]), nullptr, 0, 0, subid, pcoord, weights);

        }
        template<typename DataSet>
        int find(const V& x, double* pcoord, double* weights, const DataSet* ds, debug_info& info) {
            int subid;
            if (info.verbose) {
                info.os << "find_cell_without_locator:find with x=[" << x[0] << "," << x[1] << "," << x[2] << "], ds=" << (long)ds << "\n";
            }
            return const_cast<DataSet*>(ds)->FindCell(const_cast<S*>(&x[0]), nullptr, 0, 0, subid, pcoord, weights);
        }
    };

    template<typename V, typename S>
    struct find_cell_without_locator<V, S, typename std::enable_if<!is_double_valued<S>::value>::type> {
        template<typename DataSet>
        int find(const V& x, double* pcoord, double* weights, const DataSet* ds) {
            int subid;
            double y[3] = {x[0], x[1], x[2]};
            return const_cast<DataSet*>(ds)->FindCell(y, nullptr, 0, 0, subid, pcoord, weights);
        }
    };

    template<typename V, typename S=typename V::value_type, typename Enable=void>
    struct find_cell_with_locator {};

    template<typename V, typename S>
    struct find_cell_with_locator<V, S, typename std::enable_if<is_double_valued<S>::value>::type> {
        template<typename Locator>
        int find(const V& x, double* weights, const Locator& loc, vtkGenericCell* acell) {
            double pcoord[3];
            return loc->FindCell(const_cast<S*>(&x[0]), 0, acell, pcoord, weights);
        }
    };

    template<typename V, typename S>
    struct find_cell_with_locator<V, S, typename std::enable_if<!is_double_valued<S>::value>::type> {
        template<typename Locator>
        int find(const V& x, double* weights, const Locator& loc, vtkGenericCell* acell) {
            double pcoord[3], xcopy[3];
            std::copy(&x[0], &x[3], &xcopy[0]);
            return loc->FindCell(xcopy, 0, acell, pcoord, weights);
        }
    };

    template<typename vec>
    inline std::string to_str(const vec& p) {
        std::ostringstream os;
        if (xavier::data_traits<vec>::size() == 3)
            os << "[" << p(0) << ", " << p(1) << ", " << p(2) << "]";
        else
            os << "[" << p(0) << ", " << p(1) << "]";
        return os.str();
    }

} // anonymous namespace

template <unsigned int N>
struct locator_traits {};

template<>
struct locator_traits<3> {
    typedef vtkCellTreeLocator locator_type;
};

template<>
struct locator_traits<2> {
    typedef vtkCellLocator locator_type;
};

namespace vtk_utils {
    template< typename DataSet,
              typename T=double,
              unsigned int N=3,
              typename Vector=Eigen::Matrix<T, 3, 1>,
              typename Matrix=Eigen::Matrix<T, 3, 3>,
              typename Tuple=Eigen::Matrix<T, Eigen::Dynamic, 1>>
    class interpolator_base {
    public:
        typedef DataSet dataset_type;
        typedef T       value_type;
        typedef Vector  vector_type;
        typedef Matrix  matrix_type;
        static constexpr unsigned int dimension = N;

        interpolator_base(dataset_type* data) : m_dataset(data) {}
        interpolator_base(VTK_SMART(dataset_type) data) : m_dataset(data) {}
        const dataset_type* get_dataset() const { return m_dataset; }

        VTK_SMART(dataset_type) m_dataset;
        static inline vector_type to_vec(value_type* ptr) {
            return vector_type(*ptr, *(ptr+1), *(ptr+2));
        }
        static inline matrix_type to_mat(value_type* ptr) {
            matrix_type m;
            memcpy(&m(0,0), ptr, 9*sizeof(T));
            return m;
        }
        static inline value_type& zero(value_type& s) {
            s=0;
            return s;
        }
        static inline vector_type& zero(vector_type& v) {
            v[0] = v[1] = v[2] = 0;
            return v;
        }
        static inline matrix_type& zero(matrix_type& m) {
            std::fill(&m(0,0), &m(0,0)+9, 0);
            return m;
        }
    };

    template<typename DataSet,
             typename T=double,
             unsigned int N=3,
             typename Vector=Eigen::Matrix<T, 3, 1>,
             typename Matrix=Eigen::Matrix<T, 3, 3>,
             typename Tuple=Eigen::Matrix<T, Eigen::Dynamic, 1>,
             typename Enable=void>
    class interpolator {};

    template<typename DataSet,
             typename T,
             unsigned int N,
             typename Vector,
             typename Matrix,
             typename Tuple>
    class interpolator<DataSet, T, N, Vector, Matrix, Tuple,
                 typename std::enable_if<!needs_cell_locator<DataSet>::value>::type>
    : public interpolator_base< DataSet, T, N, Vector, Matrix, Tuple> {
    public:
        typedef interpolator_base< DataSet, T, N, Vector, Matrix, Tuple > base_type;
        typedef DataSet dataset_type;
        typedef T value_type;
        typedef Vector vector_type;
        typedef Matrix matrix_type;
        typedef Tuple  tuple_type;
        static constexpr unsigned int dimension = base_type::dimension;

        static constexpr unsigned int n_points_per_cell = (dimension == 3 ? 8 : (dimension == 2 ? 4 : 2));

        interpolator(dataset_type* data) : base_type(data), m_converter(data) {}
        interpolator(VTK_SMART(dataset_type) data) : base_type(data), m_converter(data) {}

        bool interpolation_info(const vector_type& x,
            std::vector<double>& weights, std::vector<int>& ptids,
            int& cellid) const {
                weights.resize(8);
                ptids.resize(8);
                if (!coeff(x, &weights[0], cellid)) return false;
                ptids[0] = m_converter.cellid_2_pointid(cellid);
                for (size_t i=1; i<n_points_per_cell; ++i) {
                    ptids[i] = ptids[0] + m_converter.delta(i);
                }
                return true;
        }

        bool interpolate(value_type& value, const vector_type& x) const { // scalar interpolation
            double weights[8];
            int cellid;
            if (!coeff(x, weights, cellid)) return false;
            int ptid = m_converter.cellid_2_pointid(cellid);
            vtkDataArray* scalars = this->m_dataset->GetPointData()->GetScalars();
            base_type::zero(value);
            double s;
            for (size_t i=0; i<n_points_per_cell; ++i) {
                scalars->GetTuple(ptid + m_converter.delta(i), &s);
                value += weights[i] * s;
            }
            return true;
        }

        bool interpolate(vector_type& value, const vector_type& x) const { // vector interpolation
            double weights[8];
            int cellid;

            if (!coeff(x, weights, cellid)) return false;
            int ptid = m_converter.cellid_2_pointid(cellid);
            vtkDataArray* vectors = this->m_dataset->GetPointData()->GetVectors();
            base_type::zero(value);
            double v[3];
            for (size_t i=0; i<n_points_per_cell; ++i) {
                vectors->GetTuple(ptid + m_converter.delta(i), v);
                value += weights[i] * base_type::to_vec(v);
            }
            return true;
        }

        bool interpolate(vector_type& value, const vector_type& x, bool verbose) const { // vector interpolation
            double weights[8];
            int cellid;

            debug_info dbginf;
            dbginf.verbose = verbose;

            if (verbose) {
                dbginf.os << "interpolate() called at x=[" << x[0] << "," << x[1] << "," << x[2] << "]\n";
            }

            if (!coeff(x, weights, cellid, dbginf)) {
                if (verbose) {
                    dbginf.os << "coeff returned false" << std::endl;
                    std::cout << dbginf.os.str();
                }
                return false;
            }
            int ptid = m_converter.cellid_2_pointid(cellid);

            if (verbose) {
                dbginf.os << "ptid(" << cellid << ")=" << ptid << '\n';
                dbginf.os << "coeff returned celldid=" << cellid << " and weights=[";
                std::copy(weights, weights+n_points_per_cell, std::ostream_iterator<double>(dbginf.os, ","));
                dbginf.os << "]\n";
            }
            vtkDataArray* vectors = this->m_dataset->GetPointData()->GetVectors();

            base_type::zero(value);
            if (verbose) {
                dbginf.os << "value=[" << value[0] << "," << value[1] << "," << value[2];
                dbginf.os << "]\n";
            }
            double v[n_points_per_cell];
            for (size_t i=0; i<n_points_per_cell; ++i) {
                vectors->GetTuple(ptid + m_converter.delta(i), v);
                value += weights[i] * base_type::to_vec(v);
                if (verbose) {
                    vector_type w = base_type::to_vec(v);
                    dbginf.os << "vector value #" << i << " (" << ptid+m_converter.delta(i) << ") is ["
                        << w[0] << "," << w[1] << "," << w[2];
                    dbginf.os << "]\n";
                }
            }
            if (verbose) {
                dbginf.os << "computed value=[" << value[0] << ","
                    << value[1] << ", " << value[2] << "]" << std::endl;
                std::cout << dbginf.os.str();
            }
            return true;
        }

        bool interpolate(matrix_type& value, const vector_type& x) const { // tensor interpolation
            double weights[8];
            int cellid;
            if (!coeff(x, weights, cellid)) return false;
            int ptid = m_converter.cellid_2_pointid(cellid);
            vtkDataArray* tensors = this->m_dataset->GetPointData()->GetTensors();
            double m[9];
            base_type::zero(value);
            for (size_t i=0; i<n_points_per_cell; ++i) {
                tensors->GetTuple(ptid+m_converter.delta(i), m);
                value += weights[i] * base_type::to_mat(m);
            }
            return true;
        }

        bool interpolate(tuple_type& value, const vector_type& x,
                         const std::string& field_name) const {
             double weights[8];
             int cellid;
             if (!coeff(x, weights, cellid)) return false;
             int ptid = m_converter.cellid_2_pointid(cellid);
             vtkDataArray* tuples = this->m_dataset->GetPointData()->GetAbstractArray(field_name.c_str());
             double m[100];
             int n = tuples->GetNumberOfComponents();
             value = tuple_type(n);
             for (int i=0; i<n; ++i) value(i) = 0.;

             for (size_t i=0; i<n_points_per_cell; ++i) {
                 tuples->GetTuple(ptid+m_converter.delta(i), m);
                 for (int k=0; k<n; ++k)
                    value(k) += weights[i] * m[k];
             }
             return true;
        }

    private:

        index_converter<dataset_type> m_converter;

        bool coeff(const vector_type& x, double* weights, int& cellid) const {
            double pcoord[dimension];

            find_cell_without_locator<vector_type, value_type> fc;
            cellid = fc.find(x, pcoord, weights, base_type::get_dataset());
            if (cellid < 0) {
                return false;
            }
            return true;
        }

        bool coeff(const vector_type& x, double* weights, int& cellid, debug_info& info) const {
            double pcoord[dimension];

            if (info.verbose) {
                info.os << "Coeff was called with x=[" << x[0] << "," << x[1];
                if (dimension == 3) info.os << "," << x[2];
                info.os << "]\n";
            }

            find_cell_without_locator<vector_type, value_type> fc;
            cellid = fc.find(x, pcoord, weights, base_type::get_dataset(), info);
            if (cellid < 0) {
                return false;
            }
            return true;
        }
    };

    template<typename DataSet,
             typename T,
             unsigned int N,
             typename Vector,
             typename Matrix,
             typename Tuple>
    class interpolator<DataSet, T, N, Vector, Matrix, Tuple,
        typename std::enable_if<needs_cell_locator<DataSet>::value>::type>
    : public interpolator_base<DataSet, T, N, Vector, Matrix, Tuple> {
    public:
        typedef T value_type;
        typedef interpolator_base<DataSet, T, N, Vector, Matrix, Tuple> base_type;
        typedef DataSet dataset_type;
        typedef Vector vector_type;
        typedef Matrix matrix_type;
        typedef typename locator_traits<N>::locator_type locator_type;
        static constexpr unsigned int dimension = base_type::dimension;

        interpolator(DataSet* data, bool verbose=false) : base_type(data), m_verbose(verbose) {
            m_locator = VTK_SMART(locator_type)::New();
            m_locator->SetDataSet(this->m_dataset);
            m_locator->AutomaticOn();
            m_locator->BuildLocator();
        }
        interpolator(VTK_SMART(DataSet) data, bool verbose=false) : base_type(data), m_verbose(verbose) {
            m_locator = VTK_SMART(locator_type)::New();
            m_locator->SetDataSet(this->m_dataset);
            m_locator->AutomaticOn();
            m_locator->BuildLocator();
        }

        bool interpolate(value_type& value, const vector_type& x) const {
            vtkDataArray* scalars = this->m_dataset->GetPointData()->GetScalars();
            VTK_CREATE(vtkGenericCell, acell);
            double weights[8];
            if (!coeff(x, weights, acell)) return false;

            if (m_verbose) {
                std::cerr << "point " << to_str(x) << " is located in cell ";
                acell->PrintSelf(std::cerr, vtkIndent(0));
                std::cerr << '\n';
            }

            value = 0;
            double s;
            for (size_t i=0; i<acell->GetNumberOfPoints(); ++i) {
                scalars->GetTuple(acell->GetPointId(i), &s);
                value += weights[i] * s;
            }
            return true;
        }

        bool interpolate(vector_type& value, const vector_type& x) const {
            vtkDataArray* vectors = this->m_dataset->GetPointData()->GetVectors();
            VTK_CREATE(vtkGenericCell, acell);
            double weights[8];
            if (!coeff(x, weights, acell)) return false;

            if (m_verbose) {
                std::cerr << "point " << to_str(x) << " is located in cell ";
                acell->PrintSelf(std::cerr, vtkIndent(0));
                std::cerr << '\n';
            }

            base_type::zero(value);
            double v[3];
            for (size_t i=0; i<acell->GetNumberOfPoints(); ++i) {
                vectors->GetTuple(acell->GetPointId(i), v);
                value += weights[i] * base_type::to_vec(v);
            }
            return true;
        }

        bool interpolate(matrix_type& value, const vector_type& x) const {
            vtkDataArray* tensors = this->m_dataset->GetPointData()->GetTensors();
            VTK_CREATE(vtkGenericCell, acell);
            double weights[8];
            if (!coeff(x, weights, acell)) return false;
            base_type::zero(value);
            double m[9];
            for (size_t i=0; i<acell->GetNumberOfPoints(); ++i) {
                tensors->GetTuple(acell->GetPointId(i), m);
                value += weights[i] * base_type::to_mat(m);
            }
            return true;
        }

        // thread safe versions
        bool interpolate(value_type& value, const vector_type& x, vtkGenericCell* acell) const {
            vtkDataArray* scalars = this->m_dataset->GetPointData()->GetScalars();
            double weights[8];
            if (!coeff(x, weights, acell)) return false;
            value = 0;
            double s;
            for (size_t i=0; i<acell->GetNumberOfPoints(); ++i) {
                scalars->GetTuple(acell->GetPointId(i), &s);
                value += weights[i] * s;
            }
            return true;
        }

        bool interpolate(vector_type& value, const vector_type& x, vtkGenericCell* acell) const {
            vtkDataArray* vectors = this->m_dataset->GetPointData()->GetVectors();
            double weights[8];
            if (!coeff(x, weights, acell)) return false;

            if (m_verbose) {
                std::cerr << "point " << to_str(x) << " is located in cell ";
                acell->PrintSelf(std::cerr, vtkIndent(0));
                std::cerr << '\n';
            }

            base_type::zero(value);
            double v[3];
            for (size_t i=0; i<acell->GetNumberOfPoints(); ++i) {
                vectors->GetTuple(acell->GetPointId(i), v);
                value += weights[i] * base_type::to_vec(v);
            }
            return true;
        }

        bool interpolate(matrix_type& value, const vector_type& x, vtkGenericCell* acell) const {
            vtkDataArray* tensors = this->m_dataset->GetPointData()->GetTensors();
            double weights[8];
            if (!coeff(x, weights, acell)) return false;
            base_type::zero(value);
            double m[9];
            for (size_t i=0; i<acell->GetNumberOfPoints(); ++i) {
                tensors->GetTuple(acell->GetPointId(i), m);
                value += weights[i] * base_type::to_mat(m);
            }
            return true;
        }

    private:
        VTK_SMART(locator_type) m_locator;
        bool coeff(const vector_type& x, double* weights, vtkGenericCell* acell) const {
            find_cell_with_locator<vector_type, value_type> fc;
            return fc.find(x, weights, m_locator, acell) >= 0;
        }

        bool m_verbose;
    };

    template<typename DataSet,
             typename T=double,
             unsigned int N=3,
             typename Vector=Eigen::Matrix<T, 3, 1>,
             typename Index=int,
             typename Enable=void>
    class point_locator {};

    /*
    template<typename DataSet,
             typename T,
             unsigned int N,
             typename Vector,
             typename Index>
    class point_locator<DataSet, T, N, Vector, Index,
             typename std::enable_if<!needs_cell_locator<DataSet>::value>::type> {
    public:
        typedef DataSet dataset_type;
        typedef T value_type;
        typedef Vector vector_type;
        typedef Vector point_type;
        typedef Index index_type;
        static constexpr unsigned int dimension = N;

        point_locator(dataset_type* data) : m_dataset(data), m_converter(data) {}
        point_locator(VTK_SMART(dataset_type) data) : m_dataset(data), m_converter(data) {}

        bool find_cell(const vector_type& x, std::vector<double>& weights, std::vector<index_type>& ptids) const {
            double pcoord[3];

            find_cell_without_locator<vector_type, value_type> fc;
            weights.resize(8);
            int cellid = fc.find(x, pcoord, (double*)&weights[0], m_dataset);
            if (cellid < 0)  return false;
            m_converter.cellid_2_vertids(cellid, ptids);
            return true;
        }

        VTK_SMART(dataset_type) get_dataset() const { return m_dataset; }

    private:
        VTK_SMART(dataset_type) m_dataset;
        index_converter<dataset_type> m_converter;
    };
    */

    template<typename DataSet,
             typename T,
             unsigned int N,
             typename Vector,
             typename Index>
    class point_locator<DataSet, T, N, Vector, Index,
             typename std::enable_if<needs_cell_locator<DataSet>::value>::type> {
    public:
        typedef DataSet dataset_type;
        typedef T value_type;
        typedef Vector vector_type;
        typedef Vector point_type;
        typedef Index index_type;
        typedef typename locator_traits<N>::locator_type locator_type;

        static constexpr unsigned int dimension = N;

        point_locator(dataset_type* data, bool verbose=false, bool fast=false)
            : m_dataset(data), m_verbose(verbose), m_fast(fast) {
            m_locator = VTK_SMART(locator_type)::New();
            if (verbose) m_locator->DebugOn();
            m_locator->SetDataSet(this->m_dataset);
            m_locator->AutomaticOn();
            m_locator->SetCacheCellBounds(1);
            m_locator->BuildLocator();
        }
        point_locator(VTK_SMART(dataset_type) data, bool verbose=false, bool fast=false)
            : m_dataset(data), m_verbose(verbose), m_fast(fast) {
            m_locator = VTK_SMART(typename locator_traits<dimension>::locator_type)::New();
            if (verbose) m_locator->DebugOn();
            m_locator->SetDataSet(this->m_dataset);
            m_locator->AutomaticOn();
            m_locator->BuildLocator();
        }

        bool find_cell(const vector_type& x, std::vector<double>& weights, std::vector<index_type>& ptids) const {
            VTK_CREATE(vtkGenericCell, acell);
            return find_cell(x, weights, ptids, acell);
        }

        bool find_cell(const vector_type& x, std::vector<double>& weights,
                       std::vector<index_type>& ptids, VTK_SMART(vtkGenericCell) acell) const {
            double pcoord[3];
            find_cell_with_locator<vector_type, value_type> fc;
            weights.resize(8);
            if (fc.find(x, (double*)&weights[0], m_locator, acell) < 0) return false;

            if (m_verbose) {
                std::cerr << "point " << to_str(x) << " is located in cell ";
                acell->PrintSelf(std::cerr, vtkIndent(0));
                std::cerr << '\n';
            }

            ptids.resize(acell->GetNumberOfPoints());
            weights.resize(acell->GetNumberOfPoints());
            for (int i=0; i<ptids.size(); i++) {
                ptids[i] = acell->GetPointId(i);
            }
            return true;

        }

        VTK_SMART(dataset_type) get_dataset() const { return m_dataset; }

    private:
        VTK_SMART(dataset_type) m_dataset;
        VTK_SMART(locator_type) m_locator;
        bool m_verbose;
        bool m_fast;
    };
}

#endif
