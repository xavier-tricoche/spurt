#ifndef __newAMR_HPP__
#define __newAMR_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt {
template<typename T, typename PD, typename CD>
struct AMR_quad {
    typedef fixed_vector<T, 2>        pos_type;
    typedef nvis::bounding_box<pos_type>    bounds_type;
    typedef PD                              point_data_type;
    typedef CD                              cell_data_type;
    
    bounds_type bounds() const {
        return bounds_type(pos[0], pos[2]);
    }
    
    pos_type            pos[4];
    point_data_type     point_val[4];
    cell_data_type      cell_val;
    size_t              depth;
};

struct Lt_index {
    static bool operator()(const ivec3& i0, const ivec3& i1) {
        // check depth first
        if (i0[2] < i1[2]) {
            return true;
        } else if (i0[2] > i1[2]) {
            return false;
        }
        // lexicographical order
        if (i0[0] < i1[0]) {
            return true;
        } else if (i0[0] > i1[0]) {
            return false;
        }
        return (i0[1] < i1[1]);
    }
};

// friendly wrapper around std::map<Key, Value>
template<typename D>
class AMR_data_container {
public:
    typedef ivec3                                 index_type;
    typedef D                                           data_type;
    typedef std::map<index_type, data_type, Lt_index>   map_type;
    typedef map_type::key_type                          key_type;
    typedef map_type::iterator                          iterator_type;
    typedef map_type::const_iterator                    const_iterator_type;
    
    AMR_data_container() {}
    
    bool empty() const {
        return __data.empty();
    }
    
    size_t size() const {
        return __data.size();
    }
    
    bool exists(const index_type& id) const {
        return (__data.find(id) != __data.end());
    }
    
    const D& get(const index_type& id) const {
        map_type::const_iterator __iter = __data.find(id);
        if (__iter == __data.end()) {
            throw std::runtime_error("AMR_data_container: data does not exist.");
        } else {
            return __iter->second;
        }
    }
    
    D& get(const index_type& id) {
        map_type::iterator __iter = __data.find(id);
        if (__iter == __data.end()) {
            throw std::runtime_error("AMR_data_container: key does not exist.");
        } else {
            return __iter->second;
        }
    }
    
    bool safe_get(data_type data, const index_type& id) {
        try {
            data = get(id);
        } catch(...) {
            return false;
        }
    }
    
    void add(const index_type& id, const data_type& data) {
        __data.insert(key_type(id, data));
    }
    
    void remove(const index_type& id) {
        map_type::iterator __iter = __data.find(id);
        if (__iter == __data.end()) {
            std::cerr << "AMR_data_container: Warning: key does not exist.";
        }
    }
    
    iterator_type begin() {
        return __data.begin();
    }
    const_iterator_type begin() const {
        return __data.begin();
    }
    iterator_type end() {
        return __data.end();
    }
    const_iterator_type end() const {
        return __data.end();
    }
    
private:
    map_type    __data;
};

template<typename T, typename PD, typename CD>
class AMR_root;

template<typename T, typename PD, typename CD>
class AMR_node {
public:
    typedef T                                   value_type;
    typedef AMR_node<value_type, PD, CD>        self_type;
    typedef AMR_root<T, PD, CD>                 root_type;
    typedef PD                                  point_data_type;
    typedef CD                                  cell_data_type;
    typedef fixed_vector<value_type, 2>   pos_type;
    typedef nvis::bounding_box<pos_type>        bounds_type;
    typedef self_type*                          pointer_type;
    typedef const self_type*                    const_pointer_type;
    typedef ivec3                         index_type;
    
private:
    static const index_type __corner_shift[] = {
        index_type(0,0,0), index_type(1,0,0),
        index_type(1,1,0), index_type(0,1,0)
    };
    
public:
    AMR_node(const index_type& node_index, const AMR_node* parent, const root_type* contents)
        : __index(node_index), __root(contents), __parent(parent) {
        for (int i=0 ; i<4 ; ++i) {
            __child[i] = NULL;
        }
    }
    
    ~AMR_node() {
        collapse();
    }
    
    size_t depth() const {
        return __index[2];
    }
    
    const index_type& index() const {
        return __index;
    }
    
    bounds_type bounds() const {
        return __root->cell_bounds(__index);
    }
    
    bool is_leaf() const {
        return (__children[0] == NULL);
    }
    
    const AMR_node* child(int i) const {
        return __child[i];
    }
    
    const AMR_node* parent() const {
        return __parent;
    }
    
    void collapse() const {
        if (is_leaf()) {
            return;
        }
        for (int i=0 ; i<4 ; ++i) {
            if (__child[i] != NULL) {
                delete __child[i];
            }
        }
    }
    
    AMR_quad quad() const {
        AMR_quad q;
        for (int i=0 ; i<4 ; ++i) {
            index_type id = __index + __corner_shift[i];
            q.pos[i] = __root->position(id);
            __root->point_data().safe_get(q.point_val[i], id);
        }
        __root->cell_data().safe_get(q.cell_val, __index);
        return q;
    }
    
    template<typename F>
    void split(const F& func) {
        if (!is_leaf()) {
            return;
        }
        for (int i=0 ; i<4 ; ++i) {
            index_type id(2*__index[0], 2*__index[1], depth()+1);
            id += __corner_shift[i];
            AMR_node* c = func(*this, i);
            c->__parent = this;
            c->__index = id;
            __child[i] = c;
        }
    }
    
private:
    index_type          __index;
    AMR_node*           __child[4];
    const AMR_node*     __parent;
    const root_type*    __root;
};

template<typename T, typename PD, typename CD>
class _node_leaf_iterator {
protected:
    typedef ivec3                             index_type;
    typedef AMR_node<T, PD, CD>::const_pointer_type _Base_const_ptr;
    _Base_const_ptr _M_node;
    
    inline _node_leaf_iterator(_Base_const_ptr const __N = NULL)
        : _M_node(__N) {
        // move to first reachable leaf
        if (_M_node != NULL) {
            while (!_M_node.is_leaf()) {
                _M_node = _M_node->child[0];
            }
        }
    }
    inline _node_leaf_iterator(_Base_iterator const& __THAT)
        : _M_node(__THAT._M_node) {
        // move to first reachable leaf
        if (_M_node != NULL) {
            while (!_M_node.is_leaf()) {
                _M_node = _M_node->child[0];
            }
        }
    }
    
    inline void
    _M_increment() {
        _Base_const_ptr __p = _M_node->parent;
        while (__p) {
            const index_type& id = _M_node->index();
            int i = id[0] % 2;
            int j = id[1] % 2;
            int k = i + 2*j;
            
            if (k<3) {
                _M_node = __p->child[k+1];
                break;
            } else {
                _M_node = __p;
                __p = __p->parent;
            }
        }
        _M_node = __p;
    }
    
    inline void
    _M_decrement() {
        _Base_const_ptr __p = _M_node->parent;
        while (__p) {
            const index_type& id = _M_node->index();
            int i = id[0] % 2;
            int j = id[1] % 2;
            int k = i + 2*j;
            
            if (k<3) {
                _M_node = __p->child[k+1];
                break;
            } else {
                __p = __p->parent;
            }
        }
        _M_node = __p;
        
        if (!_M_node->_M_parent) { // clearly identify the header node
            _M_node = _M_node->_M_right;
        } else if (_M_node->_M_left) {
            _Base_const_ptr x = _M_node->_M_left;
            while (x->_M_right) {
                x = x->_M_right;
            }
            _M_node = x;
        } else {
            _Base_const_ptr __p = _M_node->_M_parent;
            while (__p && _M_node == __p->_M_left) { // see below
                _M_node = __p;
                __p = _M_node->_M_parent;
            }
            if (__p) // (__p) provide undetermined behavior on rend()++ rather
                // than a seg fault, similar to standard iterator.
            {
                _M_node = __p;
            }
        }
    }
    
    template <size_t const __K, typename _Val, typename _Acc,
             typename _Dist, typename _Cmp, typename _Alloc>
    friend class KDTree;
};

};

template<typename T, typename PD, typename CD>
class AMR_root {
    static unsigned int __pow2(unsigned int n) {
        return 1 << n;
    }
    
public:
    typedef fixed_vector<T, 2>            pos_type;
    typedef PD                                  point_data_type;
    typedef CD                                  cell_data_type;
    typedef ivec3                         index_type;
    typedef nvis::bounding_box<pos_type>        bounds_type;
    typedef AMR_node<T, PD, CD>                 node_type;
    
    AMR_root(int resx, int resy, const bounds_type& bounds)
        : __res(resx, resy), __bounds(bounds) {
        __h = __bounds.size()/pos_type(__res);
        __top_layer.reserve(resx*resy);
        for (int n=0 ; n<resx*resy ; ++n) {
            int i = n%resx;
            int j = n/resx;
            __top_layer.push_back(node_type(index_type(i, j, 0), *this));
        }
    }
    
    const bounds_type& bounds() const {
        return __bounds;
    }
    
    AMR_node* get_node(const index_type& id) {
        if (it == __existing_nodes.end()) {
            return NULL;
        } else return
        }
        
    bounds_type cell_bounds(const index_type& id) const {
        int depth = id[2];
        pos_type local_h = __h/__pow2(depth);
        pos_type min = __bounds.min() + pos_type(id[0], id[1])*local_h;
        return bounds_type(min, min+local_h);
    }
    
    pos_type position(const index_type& id) const {
        int depth = id[2];
        pos_type local_h = __h/__pow2(depth);
        return __bounds.min() + pos_type(id[0], id[1])*local_h;
    }
    
    AMR_container& point_data() {
        return __point_data;
    }
    const AMR_container& point_data() {
        return __point_data;
    }
    
    AMR_container& cell_data() {
        return __cell_data;
    }
    const AMR_container& cell_data() {
        return __cell_data;
    }
    
private:
    resolution_type                             __res;
    bounds_type                                 __bounds;
    AMR_data_container<point_data_type>         __point_data;
    AMR_data_container<cell_data_type>          __cell_data;
    std::map<index_type, AMR_node*, Lt_index>   __nodes;
};

}


#endif


