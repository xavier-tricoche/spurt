#ifndef __XAVIER_KDTREE_HPP__
#define __XAVIER_KDTREE_HPP__

#include <kdtree++/kdtree.hpp>
#include <math/types.hpp>
#include <iostream>
#include <stdexcept>
#include <misc/meta_utils.hpp>

namespace spurt {

// defines a data point associating a spatial coordinate with a value
template<typename Position, typename Data>
class data_point {
public:
    typedef Position                          pos_type;
    typedef Data                              data_type;
    typedef data_traits<Position>             ptraits_type;
    typedef typename ptraits_type::value_type value_type;
    
    data_point() {}
    data_point(const pos_type& p) : m_pos(p) {}
    data_point(const pos_type& p, const data_type& d) : m_pos(p), m_data(d) {}
    data_point(const data_point& dp) : m_pos(dp.m_pos), m_data(dp.m_data) {}
    
    const pos_type& position() const {
        return m_pos;
    }

    pos_type &position()
    {
        return m_pos;
    }

    const data_type& data() const {
        return m_data;
    }
    
    data_type& data() {
        return m_data;
    }
    
    value_type distance_to(const data_point& dp) const {
        return ptraits_type::norm(dp.m_pos - m_pos);
    }
    
    value_type operator[](size_t n) const {
        return ptraits_type::value(m_pos, n);
    }
    
private:
    pos_type   m_pos;
    data_type  m_data;
};

// kdtree-based point locator - can be used for query and insertion
template <typename Position, typename Data>
struct point_locator: public KDTree::KDTree< data_traits<Position>::nrows(), 
                                             data_point<Position, Data> >
{
    typedef Position pos_type;
    typedef data_traits<pos_type> ptraits_type;
    typedef Data value_type;
    typedef data_point<Position, Data> point_type;
    typedef KDTree::KDTree<ptraits_type::nrows(), point_type> base_type;
    typedef typename base_type::const_iterator const_iterator;
    typedef typename base_type::_Region_ region_type;

    point_locator() : base_type() {}

    point_locator(const point_locator &pl) = default;

    template <typename Iterator_>
    point_locator(const Iterator_ &begin, const Iterator_ &end)
        : base_type(begin, end) {}

    point_type find_nearest_point(const pos_type &c)
    {
        if (base_type::empty())
            throw std::runtime_error("invalid query on empty tree");

        std::pair<const_iterator, value_type> found = 
            base_type::find_nearest(point_type(c), 
                                    std::numeric_limits<value_type>::max());
        return *found.first;
    }

    void find_n_nearest_points(std::list<point_type>& nn, const pos_type& c, 
                               size_t n) {
        if (base_type::empty())
            throw std::runtime_error("invalid query on empty tree");
        auto found =
            base_type::find_n_nearest(point_type(c), n);
        std::for_each(found.begin(), found.end(), [&](auto p) { 
            nn.push_back(*p.first);
        });
    }

    bool find_nearest_within_range(point_type &p, const pos_type &c, 
                                   const value_type &dist)
    {
        if (base_type::empty())
            return false;
        std::pair<const_iterator, value_type> found = 
            base_type::find_nearest(point_type(c), dist);
        if (found.first == base_type::end())
            return false;
        else
        {
            p = *found.first;
            return true;
        }
    }

    bool find_n_nearest_within_range(std::list<point_type>& nn, 
                                     const pos_type &c, size_t n, 
                                     const value_type &dist)
    {
        if (base_type::empty())
            return false;
        std::vector< std::pair<const_iterator, value_type> > found =
            base_type::find_n_nearest(point_type(c), n, dist);
        if (found.empty() || found[0].first == base_type::end())
            return false;
        else
        {
            std::for_each(found.begin(), found.end(), [&](auto p) {
                nn.push_back(*p.first);
            });
            return true;
        }
    }

    void find_within_range(std::list<point_type> &n, 
                           const pos_type &c, const value_type &dist)
    {
        base_type::find_within_range(point_type(c), dist, 
                                     std::back_inserter(n));
    }

    void find_within_range(std::list<point_type> &n, const region_type &r)
    {
        base_type::find_within_range(r, std::back_inserter(n));
    }
};

} // namespace spurt

#endif