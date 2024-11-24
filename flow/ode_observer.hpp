#pragma once

#include <vector>
#include <iostream>
#include <sstream>

#include <math/types.hpp>
#include <misc/meta_utils.hpp>

namespace spurt
{

// Observer class for boost::odeint integration API
template<typename Vector=vec3, 
         typename Scalar=typename data_traits<Vector>::value_type,
         size_t N=data_traits<Vector>::size()> 
struct Observer
{
    typedef Scalar scalar_t;
    typedef Vector vector_t;
    typedef data_traits<Vector> vector_traits_t;

    Observer(vector_t &seed, scalar_t &t, scalar_t &d, std::vector<vector_t> &c,
             std::vector<scalar_t> &ts, bool _verbose = false,
             bool save_lines = false, bool monitor = false)
        : last_p(seed), last_t(t), distance(d), curve(c), times(ts), m_verbose(_verbose), m_save_lines(save_lines), m_monitor(monitor) 
    {
        if (m_monitor || m_save_lines)
        {
            curve.push_back(seed);
            ts.push_back(t);
        }
        if (m_verbose)
        {
            std::ostringstream os;
            os << "Observer: seeding at " << to_str((seed)) << " at time " << t << std::endl;
            std::cout << os.str();
        }
    }
    
    void operator()(const vector_t &p, scalar_t t)
    {
        distance += norm(last_p - p);
        last_p = p;
        last_t = t;
        if (m_monitor || m_save_lines)
        {
            curve.push_back(p);
            times.push_back(t);
        }
        if (m_verbose)
        {
            std::ostringstream os;
            os << "\nObserver: p=" << to_str((p)) << ", t=" << t << std::endl;
            std::cout << os.str();
        }
    }

    vector_t &last_p;
    scalar_t &last_t;
    scalar_t &distance;
    std::vector<vector_t> &curve;
    std::vector<scalar_t> &times;
    bool m_verbose;
    bool m_monitor;
    bool m_save_lines;
};

} // namespace spurt