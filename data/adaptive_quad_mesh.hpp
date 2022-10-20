#ifndef __ADAPTIVE_QUAD_MESH_HPP__
#define __ADAPTIVE_QUAD_MESH_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt {
template<typename T, typename PD, typename CD>
struct AQM_quad {
    typedef nvis::fixed_vector<T, 2>        pos_type;
    typedef nvis::bounding_box<pos_type>    bounds_type;
    typedef PD                              point_data_type;
    typedef CD                              cell_data_type;
    
    bounds_type bounds() const {
        return bounds_type(pos[0], pos[2]);
    }
    
    point_data_type interpolate(const pos_type& x) const {
        T u, v, U, V;
        u = x[0];
        v = x[1];
        U = 1-u;
        V = 1-v;
        return U*V*point_val[0] + u*V*point_val[1] + U*v*point_val[3] + u*v*point_val[2];
    }
    
    nvis::fixed_vector<point_data_type, 2>
    derivative(const pos_type& x) const {
        T u, v, U, V, udot, vdot, Udot, Vdot;
        u = x[0];
        v = x[1];
        U = 1-u;
        V = 1-v;
        pos_type size = pos[2] - pos[0];
        udot = 1/size[0];
        Udot = -udot;
        vdot = 1/size[1];
        Vdot = -vdot;
        nvis::fixed_vector<point_data_type, 2> r;
        r[0] = Udot*V*point_val[0] + udot*V*point_val[1] + Udot*v*point_val[3] + udot*v*point_val[2];
        r[1] = U*Vdot*point_val[0] + u*Vdot*point_val[1] + U*vdot*point_val[3] + u*vdot*point_val[2];
        return r;
    }
    
    pos_type            pos[4];
    point_data_type     point_val[4];
    cell_data_type      cell_val;
    size_t              depth;
};

struct Lt_index {
    bool operator()(const nvis::ivec3& i0, const nvis::ivec3& i1) const {
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
class AQM_data_container {
public:
    typedef nvis::ivec3                                         index_type;
    typedef D                                                   data_type;
    typedef typename std::map<index_type, data_type, Lt_index>  map_type;
    typedef typename map_type::iterator                         iterator_type;
    typedef typename map_type::const_iterator                   const_iterator_type;
    
    AQM_data_container() {}
    
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
        typename map_type::const_iterator __iter = __data.find(id);
        if (__iter == __data.end()) {
            std::cerr << "invalid index is " << id << std::endl;
            throw std::runtime_error("AQM_data_container: data does not exist.");
        } else {
            return __iter->second;
        }
    }
    
    D& get(const index_type& id) {
        typename map_type::iterator __iter = __data.find(id);
        if (__iter == __data.end()) {
            std::cerr << "invalid index is " << id << std::endl;
            throw std::runtime_error("AQM_data_container: key does not exist.");
        } else {
            return __iter->second;
        }
    }
    
    bool safe_get(data_type data, const index_type& id) const {
        try {
            data = get(id);
        } catch(...) {
            return false;
        }
    }
    
    void add(const index_type& id, const data_type& data) {
        __data.insert(std::pair<index_type, data_type>(id, data));
    }
    
    void remove(const index_type& id) {
        typename map_type::iterator __iter = __data.find(id);
        if (__iter == __data.end()) {
            std::cerr << "AQM_data_container: Warning: key does not exist.";
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
class AQM_root;

template<typename T, typename PD, typename CD>
class AQM_node {
public:
    typedef T                                   value_type;
    typedef AQM_node<value_type, PD, CD>        self_type;
    typedef AQM_root<T, PD, CD>                 root_type;
    typedef AQM_quad<T, PD, CD>                 quad_type;
    typedef PD                                  point_data_type;
    typedef CD                                  cell_data_type;
    typedef nvis::fixed_vector<value_type, 2>   pos_type;
    typedef nvis::bounding_box<pos_type>        bounds_type;
    typedef self_type*                          pointer_type;
    typedef const self_type*                    const_pointer_type;
    typedef nvis::ivec3                         index_type;
    
private:

    static index_type id_shift(int n) {
        assert(n>=0 && n<4);
        int i = n%2;
        int j = n/2;
        return index_type(i, j, 0);
    }
    
public:
    AQM_node(const index_type& node_index, const AQM_node* parent, const root_type* contents)
        : __index(node_index), __root(contents), __parent(parent) {
        for (int i=0 ; i<4 ; ++i) {
            __child[i] = NULL;
        }
    }
    
    ~AQM_node() {
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
        return (__child[0] == NULL);
    }
    
    const self_type* child(int i) const {
        return __child[i];
    }
    
    const self_type* parent() const {
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
    
    const point_data_type& get_point_data(int i) const {
        return __root->point_data().get(index()+id_shift(i));
    }
    
    void add_point_data(int i, const point_data_type& pd) {
        return __root->point_data().add(index() + id_shift(i), pd);
    }
    
    const cell_data_type& get_cell_data() const {
        return __root->cell_data().get(__index);
    }
    
    void add_cell_data(const cell_data_type& cd) {
        return __root->cell_data().add(__index, cd);
    }
    
    quad_type quad() const {
        quad_type q;
        for (int i=0 ; i<4 ; ++i) {
            index_type id = index() + id_shift(i);
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
            id += id_shift(i);
            AQM_node* c = func(*this, i);
            c->__parent = this;
            c->__index = id;
            __child[i] = c;
        }
    }
    
    template<typename Func>
    void for_each_leaf_const(Func& f) const {
        if (is_leaf()) {
            f(this);
        } else {
            __child[0]->for_each_leaf_const(f);
            __child[1]->for_each_leaf_const(f);
            __child[3]->for_each_leaf_const(f);
            __child[2]->for_each_leaf_const(f);
        }
    }
    
    template<typename Func>
    void for_each_leaf(Func& f) {
        if (is_leaf()) {
            f(this);
        } else {
            __child[0]->for_each_leaf(f);
            __child[1]->for_each_leaf(f);
            __child[3]->for_each_leaf(f);
            __child[2]->for_each_leaf(f);
        }
    }
    
    
private:
    index_type          __index;
    AQM_node*           __child[4];
    const AQM_node*     __parent;
    const root_type*    __root;
};

// template<typename T, typename PD, typename CD>
// class AQM_node_iterator {
//     typedef nvis::ivec3                              index_type;
//     typedef AQM_node<T, PD, CD>::const_pointer_type  _Base_const_ptr;
//     _Base_const_ptr  _M_node;
//
//     AQM_node_iterator(_Base_const_ptr const __N = NULL)
//         : _M_node(__N) {
//         // move to first reachable leaf
//         if (_M_node != NULL) {
//             while (!_M_node.is_leaf()) {
//                 _M_node = _M_node->child[0];
//             }
//         }
//     }
//     AQM_node_iterator(_Base_iterator const& __THAT)
//         : _M_node(__THAT._M_node) {
//         // move to first reachable leaf
//         if (_M_node != NULL) {
//             while (!_M_node.is_leaf()) {
//                 _M_node = _M_node->child[0];
//             }
//         }
//     }
//
//  AQM_node_iterator& operator++() {
//      _M_increment();
//      return *this;
//  }
//
//  AQM_node_iterator& operator--() {
//      _M_decrement();
//      return *this;
//  }
//
//  _Base_const_ptr operator*() {
//      return _M_node;
//  }
//
//  AQM_node<T, PD, CD>& operator->() {
//      return *_M_node;
//  }
//
//
// protected:
//     inline void
//     _M_increment()
//     {
//      _Base_const_ptr __p = _M_node->parent;
//      while (__p) {
//          const index_type& id = _M_node->index();
//          int i = id[0] % 2;
//          int j = id[1] % 2;
//          int k = i + 2*j;
//
//          if (k<3) {
//              _M_node = __p->child[k+1];
//              break;
//          }
//          else {
//              _M_node = __p;
//              __p = __p->parent;
//          }
//      }
//      _M_node = __p;
//     }
//
//     inline void
//     _M_decrement()
//     {
//      _Base_const_ptr __p = _M_node->parent;
//      while (__p) {
//          const index_type& id = _M_node->index();
//          int i = id[0] % 2;
//          int j = id[1] % 2;
//          int k = i + 2*j;
//
//          if (k<3) {
//              _M_node = __p->child[k+1];
//              break;
//          }
//          else __p = __p->parent;
//      }
//      _M_node = __p;
//
//         if (!_M_node->_M_parent) // clearly identify the header node
//         {
//             _M_node = _M_node->_M_right;
//         }
//         else if (_M_node->_M_left)
//         {
//             _Base_const_ptr x = _M_node->_M_left;
//             while (x->_M_right) x = x->_M_right;
//             _M_node = x;
//         }
//         else
//         {
//             _Base_const_ptr __p = _M_node->_M_parent;
//             while (__p && _M_node == __p->_M_left) // see below
//             {
//                 _M_node = __p;
//                 __p = _M_node->_M_parent;
//             }
//             if (__p) // (__p) provide undetermined behavior on rend()++ rather
//                 // than a seg fault, similar to standard iterator.
//                 _M_node = __p;
//         }
//     }
// };

template<typename T, typename PD, typename CD>
class AQM_root {
    static unsigned int __pow2(unsigned int n) {
        return 1 << n;
    }
    
public:
    typedef nvis::fixed_vector<T, 2>            pos_type;
    typedef PD                                  point_data_type;
    typedef CD                                  cell_data_type;
    typedef nvis::ivec3                         index_type;
    typedef nvis::bounding_box<pos_type>        bounds_type;
    typedef AQM_node<T, PD, CD>                 node_type;
    
protected:
    static index_type move_up(const index_type& id, int depth) {
        int d = id[2] - depth;
        if (d<0) {
            throw std::runtime_error("invalid index in input");
        } else {
            int f = (1 << d);
            return index_type(id[0]/f, id[1]/f, depth);
        }
    }
    
public:
    AQM_root(int resx, int resy, const bounds_type& bounds)
        : __res(resx, resy), __bounds(bounds), __point_data(), __cell_data() {
        __top_layer.resize(resx*resy);
        __h = bounds.size() / nvis::vec2(resx, resy);
        std::cerr << "created an adaptive mesh with base resolution = " << resx << " x " << resy
                  << ", and bounding box = " << bounds << std::endl;
        for (int n=0 ; n<resx*resy ; ++n) {
            int i = n%resx;
            int j = n/resx;
            __top_layer[n] = new node_type(index_type(i,j,0), NULL, this);
        }
    }
    
    const bounds_type& bounds() const {
        return __bounds;
    }
    
    node_type* get_node(const index_type& id) {
        index_type _id = move_up(id, 0);
        if (_id[0] < 0 || _id[0] >= __res[0] ||
                _id[1] < 0 || _id[1] >= __res[1]) {
            throw std::runtime_error("invalid index");
        }
        
        int n = _id[0] + _id[1]*__res[0];
        node_type* _n = __top_layer[n];
        int depth = 0;
        while (depth < id[2]) {
            if (_n->is_leaf()) {
                return NULL;
            }
            index_type _base(2*_id[0], 2*_id[1], depth+1);
            _id = move_up(id, depth+1);
            index_type dif = _id - _base;
            int k = dif[0] + 2*dif[1];
            _n = _n->child[k];
            ++depth;
        }
        
        return _n;
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
    
    void add_point_data(const index_type& id, const point_data_type& pd) {
        __point_data.add(id, pd);
    }
    
    void add_cell_data(const index_type& id, const cell_data_type& cd) {
        __cell_data.add(id, cd);
    }
    
    const point_data_type& get_point_data(const index_type& id) const {
        return __point_data.get(id);
    }
    
    const cell_data_type& get_cell_data(const index_type& id) const {
        return __cell_data.get(id);
    }
    
    AQM_data_container<PD>& point_data() {
        return __point_data;
    }
    const AQM_data_container<PD>& point_data() const {
        return __point_data;
    }
    
    AQM_data_container<CD>& cell_data() {
        return __cell_data;
    }
    const AQM_data_container<CD>& cell_data() const {
        return __cell_data;
    }
    
    template<typename Func>
    void for_each_leaf_const(Func& f) const {
        for (int i=0 ; i<__top_layer.size() ; ++i) {
            __top_layer[i]->for_each_leaf_const(f);
        }
    }
    
    template<typename Func>
    void for_each_leaf(Func& f) {
        for (int i=0 ; i<__top_layer.size() ; ++i) {
            __top_layer[i]->for_each_leaf(f);
        }
    }
    
    
private:
    pos_type                                    __h;
    nvis::ivec2                                 __res;
    bounds_type                                 __bounds;
    AQM_data_container<point_data_type>         __point_data;
    AQM_data_container<cell_data_type>          __cell_data;
    std::vector<node_type*>                     __top_layer;
};

}


#endif


