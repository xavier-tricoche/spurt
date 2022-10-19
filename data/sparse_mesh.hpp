#ifndef __SPARSE_MESH_HPP__
#define __SPARSE_MESH_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <map>

namespace xavier {

template<typename T, typename PD, typename CD>
struct __Node {
    typedef T                                           value_type;
    typedef PD                                          point_data;
    typedef CD                                          cell_data;
    typedef __Node<value_type, point_data, cell_data>   self_type;
    
    __Node(const self_type* parent == NULL) : __parent(parent) {
        if (__parent) {
            __depth = __parent->depth() + 1;
        } else {
            __depth = 0;
        }
        for (int i=0 ; i<4 ; ++i) {
            __child[i] = NULL;
        }
    }
    
    void split() {
        for (int i=0 ; i<4 ; ++i) {
            __child[i] = new __Node(this);
        }
    }
    
    bool is_leaf() const {
        return (__child[0] == NULL);
    }
    
    bool is_root() const {
        return (__parent == NULL);
    }
    
    self_type* child(int i) {
        assert(i >= 0 && i < 4);
        return __child[i];
    }
    
    const self_type* child(int i) const {
        assert(i >= 0 && i < 4);
        return __child[i];
    }
    
    int depth() const {
        return __depth;
    }
    
    point_data& pdata() {
        return __pdata;
    }
    
    const point_data& pdata() const {
        return __pdata;
    }
    
    cell_data& cdata() {
        return __cdata;
    }
    
    const cell_data& cdata() const {
        return __cdata;
    }
    
    self_type*  __parent, __child[4];
    point_data  __pdata;
    cell_data   __cdata;
    int         __depth;
};



} // xavier

#endif


