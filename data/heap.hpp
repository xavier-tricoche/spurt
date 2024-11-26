#pragma once 

#include <algorithm>
#include <vector>
#include <cassert>

namespace spurt {

template< typename T,
          typename Comp = std::less<T> >
class heap : public std::vector<T> {
public:
    typedef std::vector<T> container_type;
    typedef std::vector<T> base_type;
    typedef T value_type;
    typedef Comp comparator;
    comparator m_cmp;

    heap(const comparator& cmp = comparator()) : base_type(), m_cmp(cmp) {}

    template <typename Iterator>
    heap(Iterator begin, Iterator end, const comparator& cmp = comparator()) 
        : base_type(begin, end), m_cmp(cmp) {
        std::make_heap(base_type::begin(), base_type::end(), m_cmp);
    }

    value_type top() const
    {
        return base_type::front();
    }

    value_type pop() { 
        std::pop_heap(base_type::begin(), base_type::end(), m_cmp);
        value_type r = base_type::back();
        base_type::pop_back();
        return r;
    }

    void push(const value_type& v) {
        base_type::push_back(v);
        std::push_heap(base_type::begin(), base_type::end(), m_cmp);
    }
};

template <typename T,
          typename Comp = std::less<T>>
class bounded_heap : public heap<T, Comp> {
private:
    size_t find_min() const {
        size_t mid = base_type::size()/2;
        auto min_it = std::min_element(base_type::begin() + mid, base_type::end(), base_type::m_cmp);
        return std::distance(base_type::begin(), min_it);
    }

    static size_t parent(size_t i) {
        if (i==0) return 0;
        else return (i+1)/2 - 1;
    }

    void swim(size_t from) {
        size_t pid = parent(from);
        while (from > 0 && base_type::m_cmp((*this)[pid], (*this)[from]))
        {
            std::swap((*this)[pid], (*this)[from]);
            from = pid;
            pid = parent(from);
        }
    }

    void remove_min() {
        size_t minid = find_min();
        // swap last heap entry with minimum
        std::swap((*this)[minid], (*this)[base_type::size()-1]);
        // swim up from minid position, as needed
        base_type::pop_back();
        swim(minid);
    }

public: 
    typedef heap<T, Comp> base_type;
    typedef typename base_type::container_type container_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::comparator comparator;

    bounded_heap(size_t size, const comparator& cmp=comparator()) : m_size(size), base_type(cmp) {}
    template< typename Iterator >
    bounded_heap(Iterator begin, Iterator end, size_t size, const comparator& cmp = comparator()) 
        : base_type(begin, end, cmp), m_size(size) {
        assert(container_type::size() <= m_size);
    } 

    value_type top() const {
        return base_type::top();
    }

    value_type pop() {
        return base_type::pop();
    }

    bool full() const {
        return base_type::size() == m_size;
    }

    void push(const value_type& v, bool trim_bottom=true) {
        if (container_type::size() < m_size) {
            base_type::push(v);
        }
        else if (!trim_bottom) {
            base_type::push(v);
            pop();
        }
        else {
            base_type::push(v);
            remove_min();
        }
    }

    size_t m_size;
};

} // namespace spurt