#ifndef __XAVIER_POLYNOMIAL_HPP__
#define __XAVIER_POLYNOMIAL_HPP__

// STL
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
// Boost
#include <boost/static_assert.hpp>
// arrays and such
#include <math/types.hpp>
// spurt
#include <math/combinatorial.hpp>
#include <misc/meta_utils.hpp>

namespace {

template<typename T>
inline T exponentiation_by_squaring(T x, size_t n)
{
    T r(1);
    while (n) {
        if (n % 2) { // n is odd
            r *= x;
            n--;
        }
        x *= x;
        n = n >> 1;
    }
    return r;
}
}

namespace spurt { namespace polynomial {

template<typename T>
inline T power(T x, size_t n) {
    return std::pow((double)x, (double)n);
}

template<typename, size_t, typename>
class polynomial;

template<typename Scalar, size_t Dim, 
         typename Coordinate=spurt::small_vector<Scalar, Dim> >
class monomial {
public:
    static const size_t dim = Dim;   

    typedef size_t                             degree_type;
    typedef Scalar                             value_type;
    typedef std::pair<size_t, size_t>          term_type;
    typedef Coordinate                         point_type;
    typedef monomial<Scalar, dim, Coordinate>  self_type;
    typedef std::map<size_t, size_t>           map_type;
    typedef map_type::const_iterator           const_iterator;
    typedef map_type::iterator                 iterator;
    typedef std::array<self_type, dim>         derivative_type;

    // default constructor: creates a null monomial
    monomial() : _degree(0), _constant(0), _terms() {}
    
    // range constructor
    template<typename _Iter, 
             REQ3SUBATTR(pair,_Iter::value_type,
                         integral,_Iter::value_type::first_type,
                         integral,_Iter::value_type::second_type)>
    monomial(_Iter begin, _Iter end, const value_type& constant=1)
        : _constant(constant), _degree(0) {
        if (constant != 0) {
            for (_Iter it=begin ; it!=end ; ++it) {
                _terms[it->first] = it->second;
                _degree += it->second;
            }
        }
    }
    
    // single term constructor
    monomial(const term_type& term, const value_type& constant=1)
        : _constant(constant) {
        if (constant != 0) {
            _degree = term.second;
            _terms[term.first] = term.second;
        }
    }
    
    // two term constructor
    monomial(const term_type& term1, const term_type& term2, 
             const value_type& constant=1)
        : _constant(constant), _degree(term1.second + term2.second) {
        if (constant == 0) nullify();
        else {
            _terms[term1.first] = term1.second;
            _terms[term2.first] = term2.second;
        }
    }
    
    // vector constructor
    template<typename _Vector, REQSUBATTR(integral,_Vector::value_type)>
    monomial(const _Vector& vec, const value_type& constant=1) 
        : _constant(constant), _degree(0) {
        if (constant != 0 && vec.size() > dim) {
           std::ostringstream os;
           os << "Invalid dimension (" << vec.size() << " > dim = " << dim 
              << ") in input array.\n"
              << "Caught in monomial<_Vector>(const _Vector&, const value_type&);";
           throw std::invalid_argument(os.str().c_str());
        }
        else if (constant != 0) {
            for (size_t i=0 ; i<vec.size() ; ++i) {
                if (vec[i]) {
                    _terms[i] = vec[i];
                    _degree += vec[i];
                }
            }
        }
    }
    
    // constructor from scalar
    explicit monomial(const value_type& constant) 
        : _constant(constant), _degree(0) {}
    
    // copy constructor
    template<typename OtherScalar, typename OtherCoord>
    monomial(const monomial<OtherScalar, dim, OtherCoord>& other) 
        : _constant(other._constant), _degree(other._degree),
          _terms(other._terms.begin(), other._terms.end()) {}
    
    // data access
    const value_type& constant() const {
        return _constant;
    }
    
    const_iterator begin() const {
        return _terms.begin();
    }
    
    const_iterator end() const {
        return _terms.end();
    }

    // return degree of given dimension
    degree_type degree(size_t d) const {
        const_iterator _m_it = _terms.find(d);
        if (_m_it != _terms.end()) return _m_it->second;
        else return 0;
    }
    
    // return degree of a given term
    const degree_type& degree() const {
        return _degree;
    }

    bool is_null() const {
        return _constant == 0;
    }
    
    bool operator==(const value_type& val) const {
        return _degree == 0 && _constant == val;
    }
    
    // multiply monomials
    self_type& operator*=(const term_type& term) {
        if (!is_null()) {
            iterator it = _terms.find(term.first);
            if (it != _terms.end()) it->second += term.second;
            else _terms[it->first] = it->second;
            _degree += term.second;
        }
        return *this;        
    }
    
    self_type& operator*=(const self_type& other) {
        // special case of 0 constant
        if (other._constant == 0) nullify();
        else if (!is_null()) {
            _constant *= other->_constant;
            for (const_iterator it=other._terms.begin() ;
                 it!=other._terms.end() ; ++it) {
                iterator _m_it = _terms.find(it->first);
                if (_m_it != _terms.end()) _m_it->second += it->second;
                else _terms[it->first] = it->second;
            }
            _degree += other._degree;
        }
        return *this;
    }
    
    self_type& operator*=(const value_type& constant) {
        if (constant == 0) nullify();
        else if (!is_null()) _constant *= constant;
        return *this;
    }
    
    // evaluate monomial
    value_type operator()(const point_type& x) const {
        value_type r = _constant;
        for (const_iterator it=_terms.begin() ; it!=_terms.end() ; ++it) {
            r *= power(x[it->first], it->second);
        }
        return r;
    }
    
    // return monomial corresponding to partial derivative
    self_type derivative(size_t d) const {
        if (_degree == 0 || degree(d) == 0) return self_type();
        else {
            self_type deriv(*this);
            iterator it = deriv._terms.find(d);
            if (it->second == 1) {
                deriv._terms.erase(it);
            }
            else {
                deriv._constant *= it->second;
                --it->second;
            }
            return deriv;
        }
    }
    
    // return array of partial derivatives
    derivative_type derivative() const {
        derivative_type r;
        for (size_t i=0 ; i<dim ; ++i) {
            r[i] = derivative(i);
        }
        return r;
    }
    
private:
    void nullify() {
        _constant = 0;
        _terms.clear();
        _degree = 0;
    }
    
    // friend class polynomial<value_type, dim, void>;
    
    value_type   _constant;
    map_type     _terms;
    degree_type  _degree;
};

template<typename Scalar1, size_t Dim, 
         typename Coord1=spurt::small_vector<Scalar1, Dim>, 
         typename Scalar2, 
         typename Coord2=spurt::small_vector<Scalar2, Dim> >
inline monomial<Scalar1, Dim, Coord1>
operator*(const monomial<Scalar1, Dim, Coord1>& m1, 
          const monomial<Scalar2, Dim, Coord2>& m2) {
    return monomial<Scalar1, Dim, Coord1>(m1) *= m2;
}

template<typename Scalar, size_t Dim, 
         typename Coord=spurt::small_vector<Scalar, Dim> >
inline monomial<Scalar, Dim, Coord>
operator*(const Scalar& scalar, const monomial<Scalar, Dim, Coord>& m) {
    return monomial<Scalar, Dim, Coord>(m) *= scalar;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
inline monomial<Scalar, Dim, Coord>
operator*(const std::pair<size_t, size_t> &term, const monomial<Scalar, Dim, Coord> &m)
{
    return monomial<Scalar, Dim, Coord>(m) *= term;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
inline bool
operator<(const monomial<Scalar, Dim, Coord> &m1, const monomial<Scalar, Dim, Coord> &m2)
{
    typedef monomial<Scalar, Dim, Coord> monomial_type;
    typedef typename monomial_type::const_iterator const_iterator;
    
    if (m2.is_null()) return false;
    else if (m1.degree() < m2.degree()) return true;
    else if (m1.degree() > m2.degree()) return false;
    else if (m1.constant() < m2.constant()) return true;
    else if (m1.constant() > m2.constant()) return false;
    else {
        // same degree, same constant coefficient
        const_iterator it1 = m1.begin();
        const_iterator it2 = m2.begin();
        for (; it1!=m1.end() && it2!=m2.end() ; ++it1, ++it2) {
            if (it1->first < it2->first) return true;
            else if (it1->first > it2->first) return false;
            else if (it1->second < it2->second) return true;
            else if (it1->second > it2->second) return false;
            // same degree associated with same dimension
        }
    }
    
    // these monomials are equal
    return false;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
inline bool
operator==(const monomial<Scalar, Dim, Coord> &m1, const monomial<Scalar, Dim, Coord> &m2)
{
    typedef monomial<Scalar, Dim, Coord> monomial_type;
    typedef typename monomial_type::const_iterator const_iterator;
    
    if (m1.is_null() && m2.is_null()) return true;
    else if (m1.is_null() ^ m2.is_null()) return false;
    else if (m1.degree() != m2.degree()) return false;
    else if (m1.constant() != m2.constant()) return false;
    else {
        // same degree, same constant coefficient
        const_iterator it1 = m1.begin();
        const_iterator it2 = m2.begin();
        for (; it1!=m1.end() && it2!=m2.end() ; ++it1, ++it2) {
            if (it1->first != it2->first) return false;
            else if (it1->second != it2->second) return false;
            // same degree associated with same dimension
        }
    }
    
    // these monomials are equal
    return true;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
inline bool
proportional(const monomial<Scalar, Dim, Coord> &m1, const monomial<Scalar, Dim, Coord> &m2)
{
    typedef monomial<Scalar, Dim, Coord> monomial_type;
    typedef typename monomial_type::const_iterator const_iterator;
    
    if (m1.is_null() || m2.is_null()) return false;
    else if (m1.degree() != m2.degree()) return false;
    else {
        const_iterator it1 = m1.begin();
        const_iterator it2 = m2.begin();
        for (; it1!=m1.end() && it2!=m2.end() ; ++it1, ++it2) {
            if (it1->first != it2->first) return false;
            else if (it1->second != it2->second) return false;
            // same degree associated with same dimension
        }
    }
    
    // these monomials are equal up to a constant
    return true;
}

struct partial_monomial_order {
    template <typename Scalar, size_t Dim,
              typename Coord = spurt::small_vector<Scalar, Dim>>
    bool operator()(const monomial<Scalar, Dim, Coord> &m1,
                    const monomial<Scalar, Dim, Coord> &m2) const
    {
        typedef typename monomial<Scalar, Dim, Coord>::const_iterator const_iterator;
        if (m2.is_null()) return false;
        else if (m1.degree() < m2.degree()) return true;
        else if (m1.degree() > m2.degree()) return false;
        else {
            const_iterator it1 = m1.begin();
            const_iterator it2 = m2.begin();
            for (; it1!=m1.end() && it2!=m2.end() ; ++it1, ++it2) {
                if (it1->first < it2->first) return true;
                else if (it1->first > it2->first) return false;
                else if (it1->second < it2->second) return true;
                else if (it1->second > it2->second) return false;
                // same degree associated with same dimension
            }
        }
        return false;
    }
};

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
std::ostream &operator<<(std::ostream &os, const monomial<Scalar, Dim, Coord> &m)
{
    typedef monomial<Scalar, Dim, Coord> monomial_type;
    typedef typename monomial_type::const_iterator const_iterator;
    const std::string symbols = "xyzw";
    if (m.constant() != 1 || !m.degree()) os << m.constant();
    if (monomial_type::dim < 4) {
        for (const_iterator it=m.begin() ; it!=m.end() ; ++it) {
            if (it->second>1)
                os << symbols[it->first] << '^' << it->second;
            else if (it->second==1) os << symbols[it->first];
        }
    }
    else {
        for (const_iterator it=m.begin() ; it!=m.end() ; ++it) {
            if (it->second>1)
                os << "x[" << it->first << "]^" << it->second;
            else if (it->second==1) os << "x[" << it->first << "]";
        }
    }
    return os;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
class polynomial : private std::set<monomial<Scalar, Dim, Coord>>
{
public:
    static const size_t dim = Dim;
    
    typedef size_t                               degree_type;
    typedef Scalar                               scalar_type;
    typedef Coord                                point_type;
    typedef monomial<Scalar, Dim, Coord>         monomial_type;
    typedef std::set<monomial_type>              base_type;
    typedef typename base_type::const_iterator   const_iterator;
    typedef typename base_type::iterator         iterator;
    typedef polynomial<scalar_type, dim>         self_type;
    typedef std::array<self_type, dim>           derivative_type;
    
    // default constructor: return null polynomial
    polynomial() : _degree(0) {}
    
    // range constructor
    template<typename _Iter, REQSUBTYPE(_Iter::value_type,monomial_type)>
    polynomial(_Iter begin, _Iter end) 
        : _degree(0) {
        for (_Iter it=begin ; it!=end ; ++it) {
            if (!it->is_null()) {
                base_type::insert(*it);
                _degree = std::max(it->degree(), _degree);
            }
        }
    }
    
    // copy constructor
    polynomial(const self_type& other) 
        : _degree(other._degree), 
          base_type(other.base_type::begin(), other.base_type::end()) 
        {}
    
    // monomial promotion to polynomial
    polynomial(const monomial_type& m) 
        : _degree(m.degree()) {
        if (!m.is_null()) base_type::insert(m);
    }
    
    // data access
    const size_t degree() const {
        return _degree;
    }
    
    const_iterator begin() const {
        return base_type::begin();
    }
    
    const_iterator end() const {
        return base_type::end();
    }
    
    self_type& operator+=(const polynomial& other) {
        iterator it1;
        const_iterator it2;
        for (it1=base_type::begin(), it2=other.base_type::begin() ; 
             it1!=base_type::end() && it2!=other.base_type::end() ; ) {
            if (proportional(*it1, *it2)) {
                it1->_constant += it2->_constant;
                if (it1->_constant == 0) {
                    it1 = base_type::erase(it1);
                    ++it2;
                }
                else {
                    ++it1; ++it2;
                }
            }
            else if (*it1 < *it2) {
                ++it1;
            } 
            else {
                base_type::insert(*it2);
                ++it2;
            }
        }
        if (it2 != other.base_type::end()) {
            base_type::insert(it2, other.base_type::end());
        }
        _degree = std::max(_degree, other._degree);
        return *this;
    }
    
    self_type& operator-=(const polynomial& other) {
        iterator it1;
        const_iterator it2;
        for (it1=base_type::begin(), it2=other.base_type::begin() ; 
             it1!=base_type::end() && it2!=other.base_type::end() ; ) {
            if (proportional(*it1, *it2)) {
                it1->_constant -= it2->_constant;
                if (it1->_constant == 0) {
                    it1 = base_type::erase(it1);
                    ++it2;
                }
                else {
                    ++it1; ++it2;
                }
            }
            else if (*it1 < *it2) {
                ++it1;
            } 
            else {
                self_type cpy(*it2);
                cpy *= -1;
                base_type::insert(cpy);
                ++it2;
            }
        }
        if (it2 != other.base_type::end()) {
            for (; it2!=other.base_type::end() ; ++it2) {
                self_type cpy(*it2);
                cpy *= -1;
                base_type::insert(cpy);
            }
        }
        _degree = std::max(_degree, other._degree);
        return *this;
    }
    
    self_type& operator*=(const scalar_type& s) {
        if (s == 0) {
            base_type::clear();
            _degree = 0;
        }
        else {
            for (iterator it=base_type::begin() ; it!=base_type::end() ; ++it) {
                it->_constant *= s;
            }
        }
        return *this;
    }
    
    self_type& operator/=(const scalar_type& s) {
        if (s == 0) {
            std::string msg = "Division by 0.\n";
            msg += "caught in polynomial::operator/=(const value_type&);";
            throw std::invalid_argument(msg.c_str());
        }
        else return this->operator*=(1/s);
    }
    
    self_type& operator*=(const self_type& other) {
        iterator it1;
        const_iterator it2;
        base_type new_terms;
        for (it2=other.base_type::begin() ; it2!=other.base_type::end() ; ++it2) {
            for (it1=base_type::begin() ; it1!=base_type::end() ; ++it1) {
                new_terms.insert((*it1)*(*it2));
            }
        }
        base_type::swap(new_terms);
        new_terms.clear();
        // refactorize proportional monomials
        it1=base_type::begin();
        while (it1!=base_type::end()) {
            std::pair<iterator, iterator> equal = 
                std::equal_range(it1, base_type::end(), *it1, 
                                 partial_monomial_order());
            if (std::distance(equal.first, equal.second) > 1) {
                monomial_type m(*it1++);
                for (; it1!=equal.second ; ++it1) {
                    m._constant += it1->_constant;
                }
                if (m._constant != 0) new_terms.insert(m);
            }
            else new_terms.insert(*it1++);
        }
        base_type::swap(new_terms);
        return *this;
    }
    
    scalar_type operator()(const point_type& x) const {
        scalar_type r=0;
        for (const_iterator it=base_type::begin() ; it!=base_type::end() ; ++it) {
            r += (*it)(x);
        }
        return r;
    }
    
    self_type derivative(size_t d) const {
        self_type r;
        if (_degree == 0) return r;
        for (const_iterator it=base_type::begin() ; it!=base_type::end() ;
             ++it) {
            monomial_type deriv = it->derivative(d);
            if (!deriv.is_null()) r += deriv; // adjusts degree automatically
        }
        return r;
    }
    
    derivative_type derivative() const {
        derivative_type r;
        for (size_t i=0 ; i<dim ; ++i) {
            r[i] = derivative(i);
        }
        return r;
    }
  
private:
    degree_type  _degree;
};

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
polynomial<Scalar, Dim, Coord> operator+(const polynomial<Scalar, Dim, Coord> &p1,
                                         const polynomial<Scalar, Dim, Coord> &p2)
{
    return polynomial<Scalar, Dim, Coord>(p1) += (p2);
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
polynomial<Scalar, Dim, Coord> operator-(const polynomial<Scalar, Dim, Coord> &p1,
                                         const polynomial<Scalar, Dim, Coord> &p2)
{
    return polynomial<Scalar, Dim, Coord>(p1) -= (p2);
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
polynomial<Scalar, Dim, Coord> operator*(const polynomial<Scalar, Dim, Coord> &p1,
                                         const polynomial<Scalar, Dim, Coord> &p2)
{
    return polynomial<Scalar, Dim, Coord>(p1) *= (p2);
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
polynomial<Scalar, Dim, Coord> operator*(const Scalar &s,
                                         const polynomial<Scalar, Dim, Coord> &p)
{
    return polynomial<Scalar, Dim, Coord>(p) *= s;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
std::ostream &operator<<(std::ostream &os, const polynomial<Scalar, Dim, Coord> &p)
{
    typedef polynomial<Scalar, Dim, Coord> poly_type;
    typedef typename poly_type::monomial_type  mono_type;
    typedef typename poly_type::const_iterator const_iterator;
    
    if (p.begin() == p.end()) os << "0";
    else {
        const_iterator it=p.begin();
        os << *it++;
        for (; it!=p.end() ; ++it) {
            os << " + " << *it;
        }
    }
    return os;
}

template <typename Scalar, size_t Dim,
          typename Coord = spurt::small_vector<Scalar, Dim>>
class polynomial_basis
{
public:
    static const size_t dimension = Dim;
    
    typedef Scalar                                 scalar_type;
    typedef monomial<Scalar, Dim, Coord>           monomial_type;
    typedef std::array<monomial_type, dimension>   derivative_type;
    typedef spurt::small_vector<size_t, dimension> array_type;
    
    static void 
    compute_basis(std::vector<monomial_type>& basis, 
                  const size_t order) {
        basis.clear();
        basis.push_back(monomial_type(1));
        array_type exponents;
        for (size_t k=1 ; k<=order ; ++k) {
            std::fill(exponents.begin(), exponents.end(), 0);
            _partial_basis(basis, exponents, k);
        }
    }
    
    static void 
    compute_basis_derivative(std::vector<derivative_type>& deriv,
                             const std::vector<monomial_type>& basis) {
        deriv.resize(basis.size());
        for (size_t i=0 ; i<basis.size() ; ++i) {
            for (size_t d=0 ; d<dimension ; ++d) {
                deriv[i][d] = basis[i].derivative(d);
            }
        }
    }
    
    static void compute_basis_derivative(std::vector<derivative_type>& deriv,
                                         const size_t order) {
        std::vector<monomial_type> basis;
        compute_basis(basis, order);
        compute_basis_derivative(deriv, basis);
    }
    
private:
    
    static monomial_type _make_monomial(const array_type& array) {
        typedef typename monomial_type::term_type term_type;
        
        if (all(array == 0)) return monomial_type(1);
        else {
            std::vector<term_type> terms;
            for (size_t i=0 ; i<dimension ; ++i) {
                terms.push_back(term_type(i, array[i]));
            }
            return monomial_type(terms.begin(), terms.end());
        }
    }
    
    static void _partial_basis(std::vector<monomial_type>& basis,
                               const array_type& partial,
                               const size_t k, size_t d=0) {
        if (!k) basis.push_back(_make_monomial(partial));
        else {
            array_type local(partial);
            if (d==dimension-1) {
                local[d] = k;
                basis.push_back(_make_monomial(local));
            }
            else {
                for (int i=k ; i>=0 ; --i) { // x[d]^i, 0<=i<=k
                    local[d] = i;
                    _partial_basis(basis, local, k-i, d+1);
                }
            }
        }
    }
};

// Instantiation in 0D yields invalid template specialization
template <typename Scalar, typename Coord>
class polynomial_basis<Scalar, 0, Coord>
{
};

// alternative implementation for 2D and 3D cases

template <typename Scalar, size_t Dim, typename Coord, size_t _Order >
struct alt_polynomial_basis
{
};

template <typename Scalar, size_t Dim, typename Coord>
struct constant_basis
{
    typedef Coord vector_type;
    typedef Scalar value_type;
    
    static const size_t dimension = Dim;
    static const size_t order     = 0;
    
    static void eval(std::vector<value_type>& r, const vector_type& x) {
        r.push_back(1);
    }
    static void derivative(std::vector<vector_type>& r, const vector_type& x) {
        r.push_back(vector_type(0));
    }
};

// 2D specialization
template <typename Scalar, typename Coord, size_t _Order >
struct alt_polynomial_basis<Scalar, 2, Coord, _Order>
{
    typedef Coord  vector_type;
    typedef Scalar value_type;
    
    static const size_t dimension = 2;
    static const size_t order     = _Order;
    
    static void eval(std::vector<value_type>& r, const vector_type& x) {
        alt_polynomial_basis<value_type, 2, Coord, order-1>::eval(r, x);
        for (size_t i=0 ; i<=order ; ++i) {
            r.push_back(power(x[0], order-i)*power(x[1], i));
        }
    }
    static void derivative(std::vector<vector_type>& r, const vector_type& x) {
        alt_polynomial_basis<value_type, 2, Coord, order-1>::derivative(r, x);
        for (size_t i=0 ; i<=order ; ++i) {
            value_type ddx = order-i ? 
                (order-i)*power(x[0], order-i-1)*power(x[1], i) : 0;
            value_type ddy = i ?
                i*power(x[0], order-i)*power(x[1], i-1) : 0;
            r.push_back(vector_type(ddx, ddy));
        }
    }
};

template<typename Scalar, typename Coord>
struct alt_polynomial_basis<Scalar, 2, Coord, 0> : public constant_basis<Scalar, 2, Coord> {};

// 3D specialization
template<typename Scalar, 
         typename Coord, 
         size_t _Order >
struct alt_polynomial_basis<Scalar, 3, Coord, _Order> {
    typedef Coord  vector_type;
    typedef Scalar value_type;
    
    static const size_t dimension = 3;
    static const size_t order     = _Order;
    
    static void eval(std::vector<value_type>& r, const vector_type& x) {
        alt_polynomial_basis<value_type, 3, Coord, order-1>::eval(r, x);
        for (size_t i=0 ; i<=order ; ++i) {
            for (size_t j=0 ; j<=i; ++j) {
                r.push_back(power(x[0], order-i)*
                            power(x[1], i-j)*
                            power(x[2], j));
            }
        }
    }
    static void derivative(std::vector<vector_type>& r, const vector_type& x) {
        alt_polynomial_basis<value_type, 3, Coord, order-1>::derivative(r, x);
        for (size_t i=0 ; i<=order ; ++i) {
            for (size_t j=0 ; j<=i; ++j) {
                value_type ddx = order-i ? 
                    (order-i)*power(x[0], order-i-1)*power(x[1], i-j)*
                    power(x[2], j) : 0;
                value_type ddy = i-j ?
                    (i-j)*power(x[0], order-i)*power(x[1], i-j-1)*
                    power(x[2], j) : 0;
                value_type ddz = j ?
                    j*power(x[0], order-i)*power(x[1], i-j)*
                    power(x[2], j-1) : 0;
                r.push_back(vector_type(ddx, ddy, ddz));
            }
        }
    }
};
template<typename Scalar, typename Coord>
struct alt_polynomial_basis<Scalar, 3, Coord, 0> : public constant_basis<Scalar, 3, Coord> {};

// 4D specialization
template<typename Scalar, 
         typename Coord, size_t _Order >
struct alt_polynomial_basis<Scalar, 4, Coord, _Order> {
    typedef Coord  vector_type;
    typedef Scalar                         value_type;
    
    static const size_t dimension = 4;
    static const size_t order     = _Order;
    
    static void eval(std::vector<value_type>& r, const vector_type& x) {
        alt_polynomial_basis<value_type, 4, Coord, order-1>::eval(r, x);
        for (size_t i=0 ; i<=order ; ++i) {
            for (size_t j=0 ; i+j<=order; ++j) {
                for (size_t k=0 ; i+j+k<=order ; ++k) {
                    r.push_back(power(x[0], i)*
                                power(x[1], j)*
                                power(x[2], k)*
                                power(x[3], order-i-j-k));
                }
            }
        }
    }
    static void derivative(std::vector<vector_type>& r, const vector_type& x) {
        alt_polynomial_basis<value_type, 4, Coord, order-1>::derivative(r, x);
        for (size_t i=0 ; i<=order ; ++i) {
            for (size_t j=0 ; i+j<=order; ++j) {
                for (size_t k=0 ; i+j+k<=order ; ++k) {
                    value_type ddx = i ? 
                        i*power(x[0], i-1)*
                          power(x[1], j)*
                          power(x[2], k)*
                          power(x[3], order-i-j-k) :
                        0;
                    value_type ddy = j ?
                        j*power(x[0], i)*
                          power(x[1], j-1)*
                          power(x[2], k)*
                          power(x[3], order-i-j-k) :
                        0;
                    value_type ddz = k ?
                        k*power(x[0], i)*
                          power(x[1], j)*
                          power(x[2], k-1)*
                          power(x[3], order-i-j-k) :
                        0;
                    value_type ddw = order-i-j-k ?
                        (order-i-j-k)*power(x[0], i)*
                                      power(x[1], j)*
                                      power(x[2], k)*
                                      power(x[3], order-i-j-k-1) :
                        0;
                    r.push_back(vector_type(ddx, ddy, ddz, ddw));
                }
            }
        }
    }
};

template<typename Scalar, typename Coord >
struct alt_polynomial_basis<Scalar, 4, Coord, 0> : public constant_basis<Scalar, 4, Coord> {};
    
        
} // polynomial
} // spurt


#endif