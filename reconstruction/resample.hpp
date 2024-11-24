#ifndef __RECONSTRUCTION_RESAMPLE_HPP__
#define __RECONSTRUCTION_RESAMPLE_HPP__

#include <vector>
#include <math/types.hpp>
#include <math/bounding_box.hpp>
#include <misc/meta_utils.hpp>

namespace spurt { namespace reconstruction {

template<typename Locator, typename Point = typename Locator::point_type >
void which_neighbors(std::vector<size_t>& ns, const Point& x0, 
                     double r, Locator& pl, 
                     const std::vector<Point>& points) {
    typedef Locator                                 locator_type;
    typedef Point                                   point_type;
    typedef typename locator_type::const_iterator   const_iterator;
    typedef spurt::data_traits<Point>               point_traits;

    ns.clear();
    std::vector<size_t> neighs;
    pl.find_within_range(neighs, x0, r);
    for (auto id : neighs) {
        if (point_traits::norm(points[id]-x0) < r) {
            ns.push_back(id);
        }
    }
}

template<typename Locator, typename Point = typename Locator::point_type>
size_t how_many(const Point& x0, double r, Locator& pl,
                const std::vector<Point>& points)
{
    typedef Locator                               locator_type;
    typedef Point                                 point_type;
    typedef typename locator_type::const_iterator const_iterator;
    typedef spurt::data_traits<Point>             point_traits;

    std::vector<size_t> neighbors;
    pl.find_within_range(neighbors, x0, r);
    size_t n=0;
    for (auto id : neighbors)
    {
        if (point_traits::norm(points[id] - x0) < r) ++n;
    }
    return n;
}

template<typename Point>
struct LessRelDist
{
    typedef Point point_type;
    typedef spurt::data_traits<Point> point_traits;

    LessRelDist(const point_type& x) : m_x(x) {}

    bool operator()(const point_type& x0, const point_type& x1) const
    {
        return point_traits::norm(x0-m_x) < 
               point_traits::norm(x1-m_x);
    }

    point_type m_x;
};

template< typename Locator, 
          typename Point = typename Locator::point_type >
double what_radius(const Point& x0, int& N, 
                   double maxrad, Locator& pl,
                   const std::vector<Point>& points)
{
    typedef Locator locator_type;
    typedef Point point_type;
    typedef spurt::data_traits<Point> point_traits;
    typedef typename locator_type::const_iterator const_iterator;

    std::vector<size_t> neighbors;
    pl.find_n_nearests(neighbors, x0, N+1);
    if (!neighbors.empty()) 
        return point_traits::norm(x0 - points[neighbors.back()]);
    else
        return -1;
}

} // namespace reconstruction 
} // namespace spurt

#endif