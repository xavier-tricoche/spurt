#ifndef __EDGE_HPP__
#define __EDGE_HPP__

#include <exception>
#include <functional>
#include <math/fixed_vector.hpp>

namespace spurt {
template<typename T, typename Compare = std::less<T> >
class edge : public nvis::fixed_vector<T, 2> {
    void swap(T& t0, T& t1) {
        T tmp = t1;
        t1 = t0;
        t0 = tmp;
    }
    
public:
    typedef edge<T, Compare> self_type;
    
    edge() : nvis::fixed_vector<T, 2>() {}
    edge(const T& t0, const T& t1) : nvis::fixed_vector<T, 2>(t0, t1) {
        Compare Lt;
        if (Lt(t1,t0)) {
            swap((*this)[0], (*this)[1]);
        }
    }
    
    bool operator<(const self_type& e) const {
        Compare Lt;
        if (Lt((*this)[0], e[0])) {
            return true;
        } else if (Lt(e[0], (*this)[0])) {
            return false;
        }
        return Lt((*this)[1], e[1]);
    }
};


namespace datastructure {

struct Edge {
    Edge() : _i(0), _j(0) {}
    Edge( int i, int j ) {
        if ( i > j ) {
            _i = j;
            _j = i;
        } else {
            _i = i;
            _j = j;
        }
    }
    Edge( const Edge& e ) : _i(e._i), _j(e._j) {}
    
    Edge& operator=( const Edge& e ) {
        _i = e._i;
        _j = e._j;
        return *this;
    }
    
    int operator==( const Edge& e ) const {
        return ( _i == e._i && _j == e._j );
    }
    
    int _i, _j;
};

int operator<( const Edge& e1, const Edge& e2 )
{
    return ( e1._i < e2._i ||
             ( e1._i == e2._i && e1._j < e2._j ) );
}

struct Segment {
    Segment() : _e0(), _e1() {}
    Segment( const Edge& e0, const Edge& e1 ) {
        if ( e1 < e0 ) {
            _e0 = e1;
            _e1 = e0;
        } else {
            _e0 = e0;
            _e1 = e1;
        }
    }
    Segment( const Segment& s ) : _e0(s._e0), _e1(s._e1) {}
    
    const Edge& val1() const {
        return _e0;
    }
    const Edge& val2() const {
        return _e1;
    }
    
    int operator<( const Segment& s ) const {
        return ( _e0 < s._e0 ||
                 ( _e0 == s._e0 && _e1 < s._e1 ) );
    }
    
    Edge _e0, _e1;
};

};

};

#endif
