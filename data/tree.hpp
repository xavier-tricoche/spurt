#ifndef __TREE_HPP__
#define __TREE_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt {

namespace QuadTree {
struct _Node_base {
    typedef _Node_base*         _Base_ptr;
    typedef _Node_base const*   _Base_const_ptr;
    
    _Base_ptr _M_parent;
    _Base_ptr _M_child[4];
    
    _Node_base(_Base_ptr const __PARENT = NULL,
               _Base_ptr const __CHILD0 = NULL,
               _Base_ptr const __CHILD1 = NULL,
               _Base_ptr const __CHILD2 = NULL,
               _Base_ptr const __CHILD3 = NULL)
        : _M_parent(__PARENT) {
        _M_child[0] = __CHILD0;
        _M_child[1] = __CHILD1;
        _M_child[2] = __CHILD2;
        _M_child[3] = __CHILD3;
    }
};

template <typename _Val>
struct _Node : public _Node_base {
    using _Node_base::_Base_ptr;
    typedef _Node* _Link_type;
    
    _Val _M_value;
    
    _Node(_Val const& __VALUE = _Val(),
          _Base_ptr const __PARENT = NULL,
          _Base_ptr const __LEFT = NULL,
          _Base_ptr const __RIGHT = NULL)
        : _Node_base(__PARENT, __LEFT, __RIGHT), _M_value(__VALUE) {}
        
#ifdef KDTREE_DEFINE_OSTREAM_OPERATORS
        
    template <typename Char, typename Traits>
    friend
    std::basic_ostream<Char, Traits>&
    operator<<(typename std::basic_ostream<Char, Traits>& out,
               _Node_base const& node) {
        out << &node;
        out << " parent: " << node._M_parent;
        out << "; left: " << node._M_left;
        out << "; right: " << node._M_right;
        return out;
    }
    
    template <typename Char, typename Traits>
    friend
    std::basic_ostream<Char, Traits>&
    operator<<(typename std::basic_ostream<Char, Traits>& out,
               _Node<_Val> const& node) {
        out << &node;
        out << ' ' << node._M_value;
        out << "; parent: " << node._M_parent;
        out << "; left: " << node._M_left;
        out << "; right: " << node._M_right;
        return out;
    }
    
#endif
};



template<typename T, typename D>
class base_quadtree;

template<typename T, typename D>
struct __base_iterator {
    typedef base_quadtree<T, D>*    pointer_type;
    
    pointer_type __ptr;
    
    __base_iterator(pointer_type ptr) : __ptr(ptr) {}
    
    void increment() {
        if (__ptr->child(0)) {
            while (__ptr->child(0)) {
                __ptr = __ptr->child(0);
            }
        } else {
        
        }
        
    }
};

template<typename T, typename D>
class base_quadtree {
public:
    typedef T                                   value_type;
    typedef D                                   data_type;
    typedef fixed_vector<T, 2>            vector_type;
    typedef nvis::bounding_box<vector_type>     box_type;
    typedef base_quadtree<vect_type>            self_type;
    typedef self_type*                          pointer_type;
    
    base_quadtree(const box_type& bounds) : __bounds(bounds), __depth(0) {
        for (int i=0 ; i<4 ; ++i) {
            __child[i] = NULL;
        }
    }
    
    data_type& data() {
        return __data;
    }
    const data_type& data() const {
        return __data;
    }
    
    const box_type& bounds() const {
        return __bounds;
    }
    
    int depth() const {
        return __depth;
    }
    
    self_type& child(int i) {
        assert(i>=0 && i<4);
        return *__child[i];
    }
    const self_type& child(int i) const {
        assert(i>=0 && i<4);
        return *__child[i];
    }
    
    template<typename F>
    void split(const F& data_splitter) {
        const vector_type diagonal = __bounds.size();
        vector_type __center  = __bounds.center();
        vector_type __dx(diagonal[0], 0);
        vector_type __dy(0, diagonal[1]);
        vector_type __e0 = __bounds.min() + 0.5*__dx;
        vector_type __e1 = __bounds.min() + __dx + 0.5*__dy;
        vector_type __e2 = __bounds.max() - 0.5*__dx;
        vector_type __e3 = __bounds.min() + 0.5*__dy;
        __child[0] = new self_type(box_type(__bounds.min(), __center));
        __child[1] = new self_type(box_type(__e0, __e1));
        __child[2] = new self_type(box_type(__center, __bounds.max()));
        __child[3] = new self_type(box_type(__e3, __e2));
        for (int i=0 ; i<4 ; ++i) {
            __child[i].depth = depth+1;
        }
    }
    
private:
    box_type        __bounds;
    pointer_type    __child[4];
    int             __depth;
    data_type       __data;
};


}


#endif
