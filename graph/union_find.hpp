#pragma once 

#include <vector>
#include <iterator>
#include <map>

namespace spurt {

template<typename T>
struct UnionFind {
    typedef T value_type;
    std::vector<size_t> m_parent; 
    std::vector<size_t> m_size;
    size_t m_count;
    std::map<value_type, size_t> m_value_to_index;

    size_t get_id(const value_type& p) {
        auto iter = m_value_to_index.find(p);
        if (iter == m_value_to_index.end()) {
            std::ostringstream os;
            os << "invalid value " << p << " in QuickFind::get_id";
            throw std::runtime_error(os.str());
        }
        return iter->second;
    }

    /**
     * Initializes an empty union-find data structure with
     * a sequence of values.
     * Initially, each element is in its own set.
     */
    template<typename Iterator>
    UnionFind(Iterator begin, Iterator end) {
        m_count = std::distance(begin, end);
        m_parent.resize(m_count);
        m_size.resize(m_count);
        Iterator it = begin;
        for (size_t i = 0; i < m_count; i++, ++it) {
            m_parent[i] = i;
            m_value_to_index[*it] = i;
            m_size[i] = 1;
        }
    }

    size_t count() const {
        return m_count;
    }

    /**
     * Returns the canonical element of the set containing value v
     */
    size_t find(const value_type& v) const {
        size_t n = get_id(v);
        while (n != m_parent[n]) n = m_parent[n];
        return n;
    }

    /**
     * Returns true if the two elements are in the same set.
     */
    bool connected(const value_type& p, const value_type& q) const {
        return find(p) == find(q);
    }

    /**
     * Merges the set containing element p with the set
     * containing element q.
     */
    void union(const value_type& p, const value_type& q) {
        size_t rootP = find(p);   
        size_t rootQ = find(q); 
        if (rootP == rootQ) return;
        if (m_size[rootP] < m_size[rootQ]) {
            m_parent[rootP] = rootQ;
            m_size[rootQ] += m_size[rootP];
        }
        else {
            m_parent[rootQ] = rootP;
            m_size[rootP] += m_size[rootQ];
        }
        m_count--;
    }

    void connected_sets(std::vector<std::vector<value_type>>& ccs) const {
        ccs.clear();
        std::map<int, size_t> parent_to_rank;
        for (auto kv : m_value_to_index) {
            const value_type& v = kv.first();
            size_t p = find(v);
            auto iter = parent_to_rank.find(p);
            if (iter == parent_to_rank.end()) {
                ccs.push_back(std::vector<value_type>());
                ccs.back().push_back(v);
                parent_to_rank[p] = ccs.size()-1;
            }
            else {
                ccs[iter->second].push_back(v);
            }
        }
    }
};

} // namespace spurt