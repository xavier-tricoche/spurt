#ifndef __FLOW_TIME_DEPENDENT_FIELD_HPP__
#define __FLOW_TIME_DEPENDENT_FIELD_HPP__

#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>

namespace spurt {
    template<typename Field_>
    class time_dependent_field {
    public:
        typedef Field_                             field_type;
        typedef typename field_type::vector_type   vector_type;
        typedef typename field_type::scalar_type   scalar_type;
        typedef typename field_type::point_type    point_type;

        constexpr static size_t dimension = field_type::dimension;

    private:
        constexpr static scalar_type invalid_time = std::numeric_limits<scalar_type>::min();

    public:
        time_dependent_field()
            : m_steps(), m_times(), m_deltas(), m_min_time(invalid_time),
              m_max_time(invalid_time), m_samples(0) {
        }

        time_dependent_field(const std::vector< std::shared_ptr<field_type> >& steps, const std::vector<scalar_type>& times)
            : m_steps(steps), m_times(times), m_deltas(), m_samples(0),
              m_min_time(invalid_time), m_max_time(invalid_time) {
            
            // check that the time steps are sorted and sort them if not
            if (!std::is_sorted(m_times.begin(), m_times.end())) {
                // std::cout << "reordering the time steps" << std::endl;
                
                std::vector<size_t> indices(m_times.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
                    return m_times[i1] < m_times[i2];
                });
                std::vector<scalar_type> new_times;
                std::vector< std::shared_ptr<field_type> > new_steps;
                for (size_t i=0; i<indices.size(); ++i) {
                    new_times[i] = m_times[indices[i]];
                    new_steps[i] = m_steps[indices[i]];
                }
                m_times.swap(new_times);
                m_steps.swap(new_steps);
            }

            for (size_t i=0; i<m_steps.size()-1; ++i) {
                m_deltas.push_back(m_times[i+1]-m_times[i]);
            }
            m_min_time = m_times.front();
            m_max_time = m_times.back();
        }

        time_dependent_field(const time_dependent_field& other)
            : m_steps(other.m_steps), m_times(other.m_times), m_deltas(other.m_deltas), m_samples(other.m_samples),
              m_min_time(other.m_min_time), m_max_time(other.m_max_time) {
        }

        std::pair<scalar_type, scalar_type> time_range() const {
            return std::make_pair(m_min_time, m_max_time);
        }

        vector_type operator()(const point_type& p, scalar_type& t) const {
            if (t < m_min_time || t > m_max_time) throw std::runtime_error("invalid time coordinate: " + std::to_string(t));
            else if (t == m_min_time) return (*m_steps.front())(p);
            else if (t == m_max_time) return (*m_steps.back())(p);
            else if (t == m_min_time) return (*m_steps.front())(p);
            else if (t == m_max_time) return (*m_steps.back())(p);
            size_t hi = std::distance(m_times.begin(), std::upper_bound(m_times.begin(), m_times.end(), t));
            if (hi >= m_times.size() || hi==0) throw std::runtime_error("invalid time coordinate (2): " + std::to_string(t));
            vector_type vlo = (*m_steps[hi-1])(p);      // spatial interpolation
            vector_type vhi = (*m_steps[hi])(p);        // spatial interpolation
            scalar_type u = (t-m_times[hi-1])/m_deltas[hi-1]; // linear

            std::ostringstream os;

            m_samples += 2;

            // temporal interpolation
            return (1.-u)*vlo + u*vhi;
        }
        
        // boost::odeint compliant rhs interface
        void operator()(const point_type& x, vector_type& dxdt, scalar_type t) const {
            dxdt = this->operator()(x, t);
        }

        size_t nb_samples() const { return m_samples; }
        
        void reset_counter() const { m_samples = 0; }

    private:
        std::vector< scalar_type >                  m_times, m_deltas;
        std::vector< std::shared_ptr<field_type> >  m_steps;
        scalar_type                                 m_min_time, m_max_time;
        mutable size_t                              m_samples;
    };


} // xavier

#endif
