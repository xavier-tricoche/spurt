#ifndef __AMR_HPP__
#define __AMR_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt {
template<typename T, typename PD, typename CD>
struct AMR_quad {
    typedef fixed_vector<T, 2>              pos_type;
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
    bool operator()(const ivec3& i0, const ivec3& i1) {
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
    typedef ivec3                                       index_type;
    typedef D                                           data_type;
    typedef std::map<index_type, data_type, Lt_index>   map_type;
    typedef map_type::key_type                          key_type;
    typedef map_type::iterator                          iterator_type;
    typedef map_type::const_iterator                    const_iterator_type;
    
    AMR_data_container() {}
    
    bool empty() const {
        return m_data.empty();
    }
    
    size_t size() const {
        return m_data.size();
    }
    
    bool exists(const index_type& id) const {
        return (m_data.find(id) != m_data.end());
    }
    
    const data_type& get(const index_type& id) const {
        map_type::const_iterator m_iter = m_data.find(id);
        if (m_iter == m_data.end()) {
            throw std::runtime_error("AMR_data_container: data does not exist.");
        } else {
            return m_iter->second;
        }
    }
    
    data_type& get(const index_type& id) {
        map_type::iterator m_iter = m_data.find(id);
        if (m_iter == m_data.end()) {
            throw std::runtime_error("AMR_data_container: key does not exist.");
        } else {
            return m_iter->second;
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
        m_data.insert(key_type(id, data));
    }
    
    void remove(const index_type& id) {
        map_type::iterator m_iter = m_data.find(id);
        if (m_iter == m_data.end()) {
            std::cerr << "AMR_data_container: Warning: key does not exist.";
        }
    }
    
    iterator_type begin() {
        return m_data.begin();
    }
    const_iterator_type begin() const {
        return m_data.begin();
    }
    iterator_type end() {
        return m_data.end();
    }
    const_iterator_type end() const {
        return m_data.end();
    }
    
private:
    map_type    m_data;
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
    static const index_type m_corner_shift[] = {
        index_type(0,0,0), index_type(1,0,0),
        index_type(1,1,0), index_type(0,1,0)
    };
    
public:
    AMR_node(const index_type& node_index, const root_type* contents)
        : m_index(node_index), m_root(contents) {
        m_children[0] = NULL;
    }
    
    ~AMR_node() {
        for (int i=0 ; i<4 ; ++i) {
            if (m_children[i]) {
                delete m_children[i];
            }
        }
    }
    
    size_t depth() const {
        return m_index[2];
    }
    
    bounds_type bounds() const {
        return m_root->cell_bounds(m_index);
    }
    
    bool is_leaf() const {
        return (m_children[0] == NULL);
    }
    
    AMR_quad quad() const {
        AMR_quad q;
        for (int i=0 ; i<4 ; ++i) {
            index_type id = m_index + m_corner_shift[i];
            q.pos[i] = m_root->position(id);
            m_root->point_data().safe_get(q.point_val[i], id);
        }
        m_root->cell_data().safe_get(q.cell_val, m_index);
        return q;
    }
    
    void split(self_type children[4]) {
        for (int i=0 ; i<4 ; ++i) {
            m_children[i] = &(children[i]);
        }
    }
    
    template<typename F>
    void split(const F& func) {
        for (int i=0 ; i<4 ; ++i) {
            m_children[i] = func(*this, i);
        }
    }
    
private:
    index_type          m_index;
    const root_type*    m_root;
};

template<typename T, typename PD, typename CD>
class _node_leaf_iterator {
protected:
    typedef AMR_node<T, PD, CD>::const_pointer_type _Base_const_ptr;
    m_Base_const_ptr    _M_node;
    
    inline _node_leaf_iterator(_Base_const_ptr const m_N = NULL)
        : _M_node(m_N) {
        // move to first reachable leaf
        if (_M_node != NULL) {
            while (!_M_node.is_leaf()) {
                _M_node = _M_node->children[0];
            }
        }
    }
    inline _node_leaf_iterator(_Base_iterator const& m_THAT)
        : _M_node(m_THAT._M_node) {
        // move to first reachable leaf
        if (_M_node != NULL) {
            while (!_M_node.is_leaf()) {
                _M_node = _M_node->children[0];
            }
        }
    }
    
    inline void
    _M_increment() {
        if () {
            _M_node = _M_node->_M_right;
            while (_M_node->_M_left) {
                _M_node = _M_node->_M_left;
            }
        } else {
            _Base_const_ptr m_p = _M_node->_M_parent;
            while (m_p && _M_node == m_p->_M_right) {
                _M_node = m_p;
                m_p = _M_node->_M_parent;
            }
            if (m_p) // (m_p) provide undetermined behavior on end()++ rather
                // than a seg fault, similar to standard iterator.
            {
                _M_node = m_p;
            }
        }
    }
    
    inline void
    _M_decrement() {
        if (!_M_node->_M_parent) { // clearly identify the header node
            _M_node = _M_node->_M_right;
        } else if (_M_node->_M_left) {
            _Base_const_ptr x = _M_node->_M_left;
            while (x->_M_right) {
                x = x->_M_right;
            }
            _M_node = x;
        } else {
            _Base_const_ptr m_p = _M_node->_M_parent;
            while (m_p && _M_node == m_p->_M_left) { // see below
                _M_node = m_p;
                m_p = _M_node->_M_parent;
            }
            if (m_p) // (m_p) provide undetermined behavior on rend()++ rather
                // than a seg fault, similar to standard iterator.
            {
                _M_node = m_p;
            }
        }
    }
    
    template <size_t const m_K, typename _Val, typename _Acc,
             typename _Dist, typename _Cmp, typename _Alloc>
    friend class KDTree;
};

};

template<typename T, typename PD, typename CD>
class AMR_root {
    static unsigned int m_pow2(unsigned int n) {
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
        : m_res(resx, resy), m_bounds(bounds) {
        m_h = m_bounds.size()/pos_type(m_res);
        m_top_layer.reserve(resx*resy);
        for (int n=0 ; n<resx*resy ; ++n) {
            int i = n%resx;
            int j = n/resx;
            m_top_layer.push_back(node_type(index_type(i, j, 0), *this));
        }
    }
    
    const bounds_type& bounds() const {
        return m_bounds;
    }
    
    AMR_node* get_node(const index_type& id) {
        std::set<index_type>::const_iterator it = m_existing_nodes.find(id);
    }
    
    bounds_type cell_bounds(const index_type& id) const {
        int depth = id[2];
        pos_type local_h = m_h/m_pow2(depth);
        pos_type min = m_bounds.min() + pos_type(id[0], id[1])*local_h;
        return bounds_type(min, min+local_h);
    }
    
    pos_type position(const index_type& id) const {
        int depth = id[2];
        pos_type local_h = m_h/m_pow2(depth);
        return m_bounds.min() + pos_type(id[0], id[1])*local_h;
    }
    
    AMR_container& point_data() {
        return m_point_data;
    }
    const AMR_container& point_data() {
        return m_point_data;
    }
    
    AMR_container& cell_data() {
        return m_cell_data;
    }
    const AMR_container& cell_data() {
        return m_cell_data;
    }
    
private:
    resolution_type                     m_res;
    bounds_type                         m_bounds;
    AMR_data_container<point_data_type> m_point_data;
    AMR_data_container<cell_data_type>  m_cell_data;
    std::set<index_type>                m_existing_nodes;
    std::set<index_type>                m_leaves;
};

}


#endif


