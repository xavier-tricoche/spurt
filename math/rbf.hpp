#ifndef __RBF_HPP__
#define __RBF_HPP__

#include <vector>
#include <set>
#include <list>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <array>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>

#include <data/locator.hpp>
#include <misc/meta_utils.hpp>
#include "RBFbasis.hpp"

// using Eigen SVD-based LS solution
#include <Eigen/Core>
#include <Eigen/SVD>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

namespace xavier {

namespace RBF {

// Interpolation with compactly supported radial basis function
template<typename Value, typename T, int Dim, typename Func, unsigned Order=0,
         int NVals = data_traits<Value>::size() >
class CompactSupportRBFInterpolator {
public:
    static const size_t dimension          = Dim;
    static const unsigned polynomial_order = Order;

    typedef nvis::fixed_vector<T, Dim>                        point_type;
    typedef Value                                             value_type;
    typedef data_traits<value_type>                           value_traits;
    typedef nvis::fixed_vector<Value, Dim>                    derivative_type;
    typedef T                                                 scalar_type;
    typedef Func                                              function_type;
    typedef Eigen::SparseMatrix<T>                            matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  rhs_type;
    typedef Eigen::Triplet<T>                                 triplet_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>               vector_type;
    typedef xavier::point_locator<T, int, Dim>                locator_type;
    typedef typename locator_type::point_type                 source_type;
    typedef typename polynomial::polynomial_basis<T, Dim>     polynomial_basis;
    typedef typename polynomial_basis::monomial_type          monomial_type;
    typedef typename monomial_type::derivative_type           monomial_derivative_type;

    CompactSupportRBFInterpolator(
        const std::vector<point_type>& points,
        const std::vector<value_type>& data,
        const function_type& phi,
        bool verbose = false
    )
        : _points(points), _data(data), _phi(phi) {
        // set up linear system to solve
        _size = _points.size();
        size_t total_size = _size;
        if (polynomial_order) {
            polynomial_basis::compute_basis(_poly_basis, polynomial_order);
            total_size += _poly_basis.size();
            _poly_deriv_basis.resize(_poly_basis.size());
            for (size_t i=0 ; i<_poly_basis.size() ; ++i) {
                _poly_deriv_basis[i] = _poly_basis[i].derivative();
            }
        }
        _nrhs = NVals;
        _radius = _phi.radius();

        matrix_type A(total_size, total_size);
        std::vector<triplet_type> coeff;
        rhs_type rhs(total_size, _nrhs);

        if (verbose) {
            std::cout << "building sparse RBF matrix of dimension "
                      << total_size << "...\n";
        }
        nvis::timer _timer;
        size_t counter=0;

        // insert all sites into kd-tree using efficient constructor
        std::vector<source_type> all_sites(_size);
        for (size_t i=0 ; i<all_sites.size() ; ++i) {
            all_sites[i] = source_type(_points[i], i);
        }
        _locator = locator_type(all_sites.begin(), all_sites.end());

        std::list<source_type> in_cube;
        typedef typename std::list<source_type>::const_iterator neigh_iter_t;
        for (int i=0 ; i<_size ; ++i) {
            _locator.find_within_range(in_cube, _points[i], _radius);
            for (neigh_iter_t it=in_cube.begin(); it!=in_cube.end() ; ++it) {
                int j = it->data();
                if (j < i) continue;
                if (j == i) {
                    coeff.push_back(triplet_type(i, i, _phi(0)));
                    ++counter;
                }
                else {
                    scalar_type r = nvis::norm(it->coordinate()-_points[i]);
                    if ( r < _radius) {
                        scalar_type v = _phi(r);
                        coeff.push_back(triplet_type(i, j, v));
                        coeff.push_back(triplet_type(j, i, v));
                        ++counter;
                    }
                }
            }
            for (int j=0 ; j<_nrhs ; ++j) {
                rhs(i,j) = value_traits::value(_data[i], j);
            }
        }

        if (polynomial_order) {
            for (size_t i=0 ; i<_size ; ++i) {
                std::vector<scalar_type> b(_poly_basis.size());
                for (size_t j=0 ; j<_poly_basis.size() ; ++j) {
                    b[j] = _poly_basis[j](_points[i]);
                }
                for (size_t j=0 ; j<b.size() ; ++j) {
                    coeff.push_back(triplet_type(i, _size+j, b[j]));
                    coeff.push_back(triplet_type(_size+j, i, b[j]));
                    counter += 2;
                }
            }
            for (size_t i=_size; i<A.cols() ; ++i) {
                for (size_t j=0 ; j<_nrhs ; ++j) {
                    rhs(i,j) = 0;
                }
            }
        }
        if (verbose) {
            std::cout << "\nSparse RBF matrix built in "
                      << _timer.elapsed() << " seconds\n";
            std::cout << counter << " non zero entries: "
                      << 100*counter/(_points.size()*_points.size())
                      << "%\n";
        }
        A.setFromTriplets(coeff.begin(), coeff.end());
        Eigen::SimplicialCholesky<matrix_type> cholesky(A);
        if (verbose) {
        std::fstream f("/home/xmt/tmp/matrix.txt", std::ios::out);
        f << A;
        f.close();
        std::cout << "matrix exported to file\n";
            std::cout << "starting Cholesky solver...\n";
        }
        _timer.restart();
        rhs_type sol = cholesky.solve(rhs);
        std::cout << "Cholesky solver took " << _timer.elapsed()
                  << " seconds\n";
        double err = (A*sol - rhs).norm();
        std::cout << "residual = " << err << " (" << err/rhs.norm() << ")\n";
        _weights.resize(sol.rows());
        for (int i=0 ; i<sol.rows() ; ++i) {
            for (int j=0 ; j<_nrhs ; ++j) {
                value_traits::value(_weights[i], j) = sol(i,j);
            }
        }
    }

    CompactSupportRBFInterpolator(
        const std::vector<point_type>& points,
        const std::vector<value_type>& data,
        const std::vector<value_type>& weights,
        const function_type& phi,
        bool verbose = false
    )
        : _size(points.size()), _points(points), _data(data),
          _weights(weights), _phi(phi) {
        if (polynomial_order) {
            polynomial_basis::compute_basis(_poly_basis, polynomial_order);
            _poly_deriv_basis.resize(_poly_basis.size());
            for (size_t i=0 ; i<_poly_basis.size() ; ++i) {
                _poly_deriv_basis[i] = _poly_basis[i].derivative();
            }
        }
        _radius = _phi.radius();
    }

    value_type operator()(const point_type& x) const {
        typedef typename locator_type::const_iterator   const_iterator;
        typename std::list<source_type>::const_iterator it;

        value_type val;
        value_traits::assign(val, 0);
        std::list<source_type> in_cube;

        _locator.find_within_range(in_cube, x, _radius);
        for (it=in_cube.begin(); it!=in_cube.end() ; ++it) {
            scalar_type r = nvis::norm(it->coordinate()-x);
            if ( r < _radius) {
                val += _phi(r)*_weights[it->data()];
            }
        }
        if (polynomial_order) {
            for (size_t i=0 ; i<_poly_basis.size() ; ++i) {
                val += _poly_basis[i](x)*_weights[i+_size];
            }
        }

        return val;
    }

    derivative_type derivative(const point_type& x) const {
        typedef typename locator_type::const_iterator   const_iterator;
        typename std::list<source_type>::const_iterator it;

        derivative_type nabla_f(0);
        std::list<source_type> in_cube;

        _locator.find_within_range(in_cube, x, _radius);
        for (it=in_cube.begin(); it!=in_cube.end() ; ++it) {
            scalar_type r = nvis::norm(it->coordinate()-x);
            if ( r < _radius && r > 0) {
                point_type nabla_r = (x - it->coordinate())/r;
                scalar_type phi_prime = _phi.derivative(r);
                for (int j=0 ; j<dimension ; ++j) {
                    nabla_f[j] += phi_prime*nabla_r[j]*_weights[it->data()];
                }
            }
        }
        if (polynomial_order) {
            for (size_t i=0 ; i<_poly_deriv_basis.size() ; ++i) {
                for (size_t j=0 ; j<dimension ; ++j) {
                    nabla_f[j] += _poly_deriv_basis[i][j](x)*_weights[i+_size];
                }
            }
        }

        return nabla_f;
    }

    const std::vector<point_type>& points() const {
        return _points;
    }
    const std::vector<value_type>& data() const {
        return _data;
    }
    const std::vector<value_type>& weights() const {
        return _weights;
    }

private:
    size_t                                _size;
    size_t                                _nrhs;
    std::vector<value_type>               _weights;
    const std::vector<point_type>&        _points;
    const std::vector<value_type>&        _data;
    const function_type&                  _phi;
    scalar_type                           _radius;
    mutable locator_type                  _locator;
    std::vector<monomial_type>            _poly_basis;
    std::vector<monomial_derivative_type> _poly_deriv_basis;
};

// Interpolation with radial basis function with infinite support
template< typename Value, typename T, size_t Dim, typename Func,
          typename Solver =
          typename Eigen::ColPivHouseholderQR<
          typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >,
          unsigned Order=0 >
class InfiniteSupportRBFInterpolator {
public:
    typedef nvis::fixed_vector<T, Dim>                        point_type;
    typedef Value                                             value_type;
    typedef data_traits<value_type>                           value_traits;
    typedef nvis::fixed_vector<Value, Dim>                    derivative_type;
    typedef T                                                 scalar_type;
    typedef Func                                              function_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  matrix_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>               vector_type;
    typedef Solver                                            solver_type;
    typedef Eigen::JacobiSVD<matrix_type>                     svd_solver_type;
    typedef typename polynomial::polynomial_basis<T, Dim>     polynomial_basis;
    typedef typename polynomial_basis::monomial_type          monomial_type;
    typedef typename monomial_type::derivative_type           monomial_derivative_type;

    static const size_t          dimension = Dim;
    static const unsigned polynomial_order = Order;

    InfiniteSupportRBFInterpolator(
        const std::vector<point_type>& points,
        const std::vector<value_type>& data,
        const function_type& phi,
        bool verbose = false
    )
    : _points(points), _data(data), _phi(phi), _verbose(verbose) {

        // set up linear system to solve
        _size = _points.size();
        size_t total_size = _size;
        if (polynomial_order) {
            polynomial_basis::compute_basis(_poly_basis, polynomial_order);
            total_size += _poly_basis.size();
            _poly_deriv_basis.resize(_poly_basis.size());
            for (size_t i=0 ; i<_poly_basis.size() ; ++i) {
                _poly_deriv_basis[i] = _poly_basis[i].derivative();
            }
            total_size += _poly_basis.size();
        }
        _nrhs = value_traits::size();

        matrix_type A(total_size, total_size);
        matrix_type rhs(total_size, _nrhs);

        if (verbose) {
            std::cout << "Building dense RBF matrix of dimension "
            << total_size << "...\n";
        }
        nvis::timer _timer;
        // set value constraints
        for (size_t i=0 ; i<_size ; ++i) {
            // A(i,j) = \phi(||xi-xj||)
            // rhs(i) = fi
            A(i,i) = _phi(0);
            if (verbose) {
                std::cout << "\r" << i+1 << " / " << _size
                << " (" << 100.*(float)(i+1)/(float)_size
                << "%)          \r"
                << std::flush;
            }
            for (size_t j=i+1 ; j<_size ; ++j) {
                scalar_type r = nvis::norm(_points[i]-_points[j]);
                scalar_type v = _phi(r);
                A(i,j) = v;
                A(j,i) = v;
            }
            for (size_t j=0 ; j<_nrhs ; ++j) {
                rhs(i,j) = value_traits::value(_data[i], j);
            }
        }
        if (polynomial_order) {
            for (size_t i=0 ; i<_size ; ++i) {
                for (size_t j=0 ; j<_poly_basis.size() ; ++j) {
                    A(i,_size+j) = _poly_basis[j](_points[i]-_points[0]);
                    A(_size+j,i) = A(i,_size+j);
                }
            }
            for (size_t i=_size; i<A.cols() ; ++i) {
                A(i,i) = 0;
                for (size_t j=i+1 ; j<A.cols() ; ++j) {
                    A(i,j) = 0;
                    A(j,i) = 0;
                }
                for (size_t j=0 ; j<_nrhs ; ++j) {
                    rhs(i,j) = 0;
                }
            }
        }

        if (verbose) {
            std::cout << "\nDense RBF matrix built in " << _timer.elapsed()
                      << " seconds\n";
            std::fstream fs("~/tmp/matrix.txt", std::ios::out);
            fs << A;
            fs.close();
            fs.open("~/tmp/rhs.txt", std::ios::out);
            fs << rhs;
            fs.close();
            std::cout << "starting solver...\n";
        }

        _timer.restart();
        solver_type solver(A);
        matrix_type sol = solver.solve(rhs);
        if (verbose) {
            std::cout << "solver took " << _timer.elapsed() << " seconds\n";
            // check what we got
            double residual = (A*sol - rhs).norm();
            std::cout << "residual=" << residual
                      << " (" << residual/rhs.norm() << ")\n";
            std::fstream f("/home/xmt/tmp/matrix.txt", std::ios::out);
        f << A;
        f.close();
        std::cout << "matrix exported to file\n";
        }
        _weights.resize(sol.rows());
        for (size_t i=0 ; i<sol.rows() ; ++i) {
            for (size_t j=0 ; j<_nrhs ; ++j) {
                value_traits::value(_weights[i], j) = sol(i,j);
            }
        }
    }

    InfiniteSupportRBFInterpolator(
        const std::vector<point_type>& points,
        const std::vector<value_type>& data,
        const std::vector<value_type>& weights,
        const function_type& phi,
        bool verbose = false
    )
        : _points(points), _data(data), _weights(weights), _phi(phi),
          _verbose(verbose) {

        _size = _points.size();
        _nrhs = value_traits::size();
        polynomial_basis::compute_basis(_poly_basis, polynomial_order);
        for (size_t i=0 ; i<_poly_basis.size() ; ++i) {
            _poly_deriv_basis[i] = _poly_basis[i].derivative();
        }
    }

    value_type operator()(const point_type& x) const {
        value_type val(0);

        for (int i=0 ; i<_size ; ++i) {
            scalar_type r = nvis::norm(_points[i]-x);
            val += _phi(r)*_weights[i];
        }

        if (polynomial_order) {
            for (int i=0 ; i<_poly_basis.size() ; ++i) {
                val += _poly_basis[i](x-_points[0])*_weights[i+_size];
            }
        }

        return val;
    }

    derivative_type derivative(const point_type& x) const {
        derivative_type nabla_f(0);
        // \nabla f (x) = \sum_j \omega_j \nabla \phi_j (x)
        // \nabla \phi_j = \nabla r(x,xj) \phi'(r(x,xj))
        // \nabla r(x,xj) = (x-xj)/r(x,xj)
        for (int i=0 ; i<_size ; ++i) {
            scalar_type r = nvis::norm(_points[i]-x);
        if (r==0) continue;
            point_type nabla_r = (x - _points[i])/r;
            scalar_type phi_prime = _phi.derivative(r);
            for (int j=0 ; j<dimension ; ++j) {
                nabla_f[j] += phi_prime*nabla_r[j]*_weights[i];
            }
        }
        if (polynomial_order) {
            for (size_t i=0 ; i<_poly_deriv_basis.size() ; ++i) {
                for (size_t j=0 ; j<dimension ; ++j) {
                    nabla_f[j] += _poly_deriv_basis[i][j](x-_points[0])*_weights[i+_size];
                }
            }
        }

        return nabla_f;
    }

    const std::vector<point_type>& points() const {
        return _points;
    }
    const std::vector<value_type>& data() const {
        return _data;
    }
    const std::vector<value_type>& weights() const {
        return _weights;
    }

private:
    size_t                                _size;
    size_t                                _nrhs;
    std::vector<value_type>               _weights;
    const std::vector<point_type>&        _points;
    const std::vector<value_type>&        _data;
    const function_type&                  _phi;
    bool                                  _verbose;
    std::vector<monomial_type>            _poly_basis;
    std::vector<monomial_derivative_type> _poly_deriv_basis;
};

} // namespace RBF

} // namespace xavier


#endif
