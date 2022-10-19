#ifndef __FLOW_TIME_DEPENDENT_FIELD_HPP__
#define __FLOW_TIME_DEPENDENT_FIELD_HPP__

#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <memory>

#include <VTK/vtk_utils.hpp>
#include <VTK/vtk_interpolator.hpp>
#include <vtkGenericCell.h>

#include "math/fixed_vector.hpp"
#include <math/dopri5.hpp>

namespace {

    template<typename Scalar>
    int find_time(Scalar t, const std::vector<Scalar>& times) {
        if (t < times.front() || t > times.back()) return -1;
        else if (t == times.front()) {
            return 0;
        }
        else if (t == times.back()) {
            return times.size()-2;
        }
        else {
            size_t hi = std::distance(times.begin(), std::upper_bound(times.begin(), times.end(), t));
            if (hi >= times.size() || hi==0) return -1;
            else return hi-1;
        }
    }
}

namespace xavier {

    template<typename Field_>
    class arbitrary_time_dependent_field {
    public:
        typedef Field_                             field_type;
        typedef typename field_type::scalar_type   scalar_type;
        typedef typename field_type::point_type    point_type;
        typedef typename field_type::value_type    value_type;

        constexpr static size_t dimension = field_type::dimension;

    private:
        constexpr static scalar_type invalid_time = std::numeric_limits<scalar_type>::min();

    public:
        arbitrary_time_dependent_field()
            : m_steps(), m_times(), m_deltas(), m_min_time(invalid_time),
              m_max_time(invalid_time), m_samples(0) {
            // std::cout << "default constructor @" << (void*)this << std::endl;
        }

        arbitrary_time_dependent_field(const std::vector< std::shared_ptr<field_type> >& steps,
                                       const std::vector<scalar_type>& times)
            : m_steps(steps), m_times(times), m_deltas(), m_samples(0),
              m_min_time(invalid_time), m_max_time(invalid_time) {

            // std::cout << "time-dependent field constructor @" << (void*)this << std::endl;

            // check that the time steps are sorted and sort them if not
            if (!std::is_sorted(m_times.begin(), m_times.end())) {
                std::cout << "reordering the time steps" << std::endl;

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

        arbitrary_time_dependent_field(const arbitrary_time_dependent_field& other)
            : m_steps(other.m_steps), m_times(other.m_times), m_deltas(other.m_deltas),
              m_samples(other.m_samples),
              m_min_time(other.m_min_time), m_max_time(other.m_max_time) {
            // std::cout << "time-dependent field copy constructor @" << (void*)this
                      // << " from @" << (void*)&other << ", inherited m_samples: "
                      // << other.m_samples << std::endl;
            }

        std::pair<scalar_type, scalar_type> time_range() const {
            return std::make_pair(m_min_time, m_max_time);
        }

        value_type operator()(const point_type& p, scalar_type& t) const {
            int idx = find_time(t, m_times);
            if (idx == -1) throw std::runtime_error("invalid time coordinate: " + std::to_string(t));
            value_type vlo = (*m_steps[idx])(p);      // spatial interpolation
            value_type vhi = (*m_steps[idx+1])(p);    // spatial interpolation
            scalar_type u = (t-m_times[idx])/m_deltas[idx+1]; // linear
            if (true) {
                std::ostringstream os;
                os << "interpolation at " << p << " at time step #" << idx << " returned " << vlo << '\n';
                os << "interpolation at " << p << " at time step #" << idx+1 << " returned " << vhi << '\n';
                std::cout << os.str() << std::flush;
            }

            std::ostringstream os;

            // os << "\thi=" << hi << " before increment, counter=" << m_samples << ")\n";
            m_samples += 2;

            // temporal interpolation
            // os << "\tfield(" << p[0] << ", " << p[1] << ", " << p[2] << ", " << t << ")=" << (1.-u)*vlo + u*vhi << " (counter=" << m_samples << ") @=" << (void*)this << std::endl;
            // std::cout << os.str();
            return (1.-u)*vlo + u*vhi;
        }

        size_t nb_samples() const { return m_samples; }
        void reset_counter() const { /*std::cout << "resetting counter" << std::endpoint_locatorl;*/ m_samples = 0; }

    private:
        std::vector< scalar_type >                  m_times, m_deltas;
        std::vector< std::shared_ptr<field_type> >  m_steps;
        scalar_type                                 m_min_time, m_max_time;
        mutable size_t                              m_samples;
    };

    template<typename PointLocator, typename Container, typename Value=typename Container::value_type>
    class fixed_mesh_time_dependent_field {
    public:
        typedef PointLocator locator_type;
        typedef typename locator_type::point_type point_type;
        typedef typename locator_type::value_type scalar_type;
        typedef typename locator_type::index_type index_type;
        typedef typename locator_type::dataset_type dataset_type;
        typedef Container container_type;
        typedef Value value_type;

    private:
        static constexpr scalar_type invalid_time = std::numeric_limits<scalar_type>::min();

        value_type interpolate(index_type time_step,
                               const std::vector<scalar_type>& weights,
                               const std::vector<index_type>& ptids) const {
            std::shared_ptr<const container_type> values = m_steps[time_step];

            value_type v = value_type(weights.front()*(*values)[ptids.front()]);
            if (m_verbose) {
                std::cout << "ptid[0]=" << ptids[0] << ", weight[0]=" << weights.front()
                          << ", value[0]=" << to_str((*values)[ptids.front()]) << '\n';
            }
            for (int i=1; i<weights.size(); i++) {
                v += value_type(weights[i]*(*values)[ptids[i]]);
                if (m_verbose) {
                    std::cout << "ptid[" << i << "]=" << ptids[i] << ", weight["
                              << i << "]=" << weights[i] << ", value[" << i << "]="
                              << to_str((*values)[ptids[i]]) << '\n';
                }
            }
            return v;
        }

    public:
        fixed_mesh_time_dependent_field(bool verbose=false)
            : m_times(), m_steps(), m_locator(nullptr), m_verbose(verbose) {}

        fixed_mesh_time_dependent_field(std::shared_ptr<const locator_type> locator,
                const std::vector< std::shared_ptr<container_type> >& steps,
                const std::vector< scalar_type >& times,
                bool verbose = false)
                : m_locator(locator), m_steps(steps), m_times(times),
                  m_verbose(verbose) {
            // check that the time steps are sorted and sort them if not
            if (!std::is_sorted(m_times.begin(), m_times.end())) {
                if (m_verbose) {
                    std::cout << "reordering the time steps" << std::endl;
                }

                std::vector<size_t> indices(m_times.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
                    return m_times[i1] < m_times[i2];
                });
                std::vector<scalar_type> new_times;
                std::vector< std::shared_ptr<container_type> > new_steps;
                for (size_t i=0; i<indices.size(); ++i) {
                    new_times[i] = m_times[indices[i]];
                    new_steps[i] = m_steps[indices[i]];
                }
                m_times.swap(new_times);
                m_steps.swap(new_steps);
            }
        }

        fixed_mesh_time_dependent_field(const fixed_mesh_time_dependent_field& other)
                : m_locator(other.m_locator), m_steps(other.m_steps),
                m_times(other.m_times), m_verbose(other.m_verbose) {}

        std::pair<scalar_type, scalar_type> time_range() const {
            return std::pair<scalar_type, scalar_type>(m_times.front(), m_times.back());
        }

        value_type operator()(const point_type& p, scalar_type& t) const {
            VTK_CREATE(vtkGenericCell, a_cell);
            return this->operator()(p, t, a_cell);
        }

        value_type operator()(const point_type& p, scalar_type& t,
                              VTK_SMART(vtkGenericCell) a_cell) const {
            std::vector<scalar_type> weights(8);
            std::vector<index_type> ptids(8);
            if (!m_locator->find_cell(p, weights, ptids, a_cell)) {
                std::ostringstream os;
                os << "fixed_mesh_time_dependent_field::interpolate: invalid location: " << to_str(p) << ", t=" << t;
                if (m_verbose) {
                    std::cerr << os.str() << std::flush;
                }
                throw std::runtime_error(os.str());
            }

            index_type idx = find_time(t, m_times);
            if (idx == -1) {
                std::ostringstream os;
                os << "fixed_mesh_time_dependent_field::interpolate: invalid time coordinate: " << t << ", p=" << to_str(p);
                throw std::runtime_error(os.str());
            }
            else {
                value_type vlo = interpolate(idx, weights, ptids);      // spatial interpolation
                value_type vhi = interpolate(idx+1, weights, ptids);    // spatial interpolation
                scalar_type u = (t-m_times[idx])/(m_times[idx+1]-m_times[idx]); // linear time interpolation

                if (m_verbose) {
                    std::ostringstream os;
                    os << "interpolation at " << to_str(p) << " at time step #" << idx << " returned " << to_str(vlo) << '\n';
                    os << "interpolation at " << to_str(p) << " at time step #" << idx+1 << " returned " << to_str(vhi) << '\n';
                    std::cout << os.str() << std::flush;
                }
                vlo *= (1-u);
                vhi *= u;
                return value_type(vlo + vhi);
            }
        }

        std::shared_ptr<const locator_type> get_locator() const { return m_locator; }

        std::vector< std::shared_ptr<container_type> > m_steps;
        std::vector< scalar_type > m_times;
        std::shared_ptr< const locator_type > m_locator;
        bool m_verbose;
    };

    template<typename DataSet, typename T, unsigned int N,  typename Container, typename Value=typename Container::value_type>
    class structured_mesh_time_dependent_field {
    public:
        typedef DataSet                                  dataset_type;
        typedef vtk_utils::interpolator<DataSet, T, N>   interpolator_type;
        typedef typename interpolator_type::value_type   scalar_type;
        typedef typename interpolator_type::vector_type  point_type;
        typedef typename interpolator_type::vector_type  vector_type;
        typedef typename interpolator_type::matrix_type  matrix_type;
        typedef Container container_type;
        typedef Value value_type;

    private:
        static constexpr scalar_type invalid_time = std::numeric_limits<scalar_type>::min();

        bool is_out_of_bounds(int cellid) const {
            if (m_out_of_bounds.empty()) return false;
            else return m_out_of_bounds[cellid];
        }

        value_type interpolate(int time_step,
                               const std::vector<scalar_type>& weights,
                               const std::vector<int>& ptids) const
        {
            std::shared_ptr<const container_type> values = m_steps[time_step];

            value_type v = value_type(weights.front()*(*values)[ptids.front()]);
            if (m_verbose) {
                std::cout << "ptid[0]=" << ptids[0] << ", weight[0]=" << weights.front()
                          << ", value[0]=" << to_str((*values)[ptids.front()]) << '\n';
            }
            for (int i=1; i<weights.size(); i++) {
                v += value_type(weights[i]*(*values)[ptids[i]]);
                if (m_verbose) {
                    std::cout << "ptid[" << i << "]=" << ptids[i] << ", weight["
                              << i << "]=" << weights[i] << ", value[" << i << "]="
                              << to_str((*values)[ptids[i]]) << '\n';
                }
            }
            return v;
        }

    public:
        structured_mesh_time_dependent_field(bool verbose=false)
            : m_times(), m_steps(), m_verbose(verbose), m_out_of_bounds() {}

        structured_mesh_time_dependent_field(dataset_type* dataset,
                const std::vector< std::shared_ptr<container_type> >& steps,
                const std::vector< scalar_type >& times,
                const std::vector< bool >& out_of_bounds,
                bool verbose = false)
                : m_interpolator(new interpolator_type(dataset)), m_steps(steps), m_times(times),
                  m_verbose(verbose) , m_out_of_bounds(out_of_bounds) {
            // check that the time steps are sorted and sort them if not
            if (!std::is_sorted(m_times.begin(), m_times.end())) {
                if (m_verbose) {
                    std::cout << "reordering the time steps" << std::endl;
                }

                std::vector<size_t> indices(m_times.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
                    return m_times[i1] < m_times[i2];
                });
                std::vector<scalar_type> new_times;
                std::vector< std::shared_ptr<container_type> > new_steps;
                for (size_t i=0; i<indices.size(); ++i) {
                    new_times[i] = m_times[indices[i]];
                    new_steps[i] = m_steps[indices[i]];
                }
                m_times.swap(new_times);
                m_steps.swap(new_steps);
            }
        }

        std::pair<scalar_type, scalar_type> time_range() const {
            return std::pair<scalar_type, scalar_type>(m_times.front(), m_times.back());
        }

        value_type operator()(const point_type& p, scalar_type& t) const {
            std::vector<scalar_type> weights(8);
            std::vector<int> ptids(8);
            int cellid;
            if (!m_interpolator->interpolation_info(p, weights, ptids, cellid) ||
                is_out_of_bounds(cellid)) {
                std::ostringstream os;
                os << "structured_mesh_time_dependent_field::interpolate: invalid location: " << to_str(p) << ", t=" << t;
                if (m_verbose) {
                    std::cerr << os.str() << std::flush;
                }
                throw std::runtime_error(os.str());
            }

            int idx = find_time(t, m_times);
            if (idx == -1) {
                std::ostringstream os;
                os << "structured_mesh_time_dependent_field::interpolate: invalid time coordinate: " << t << ", p=" << to_str(p);
                throw std::runtime_error(os.str());
            }
            else {
                value_type vlo = interpolate(idx, weights, ptids);      // spatial interpolation
                value_type vhi = interpolate(idx+1, weights, ptids);    // spatial interpolation
                scalar_type u = (t-m_times[idx])/(m_times[idx+1]-m_times[idx]); // linear time interpolation

                if (m_verbose) {
                    std::ostringstream os;
                    os << "interpolation at " << to_str(p) << " at time step #" << idx << " returned " << to_str(vlo) << '\n';
                    os << "interpolation at " << to_str(p) << " at time step #" << idx+1 << " returned " << to_str(vhi) << '\n';
                    std::cout << os.str() << std::flush;
                }
                vlo *= (1-u);
                vhi *= u;
                return value_type(vlo + vhi);
            }
        }

        std::shared_ptr<interpolator_type> get_interpolator() {
            return m_interpolator;
        }

        std::shared_ptr< interpolator_type > m_interpolator;
        std::vector< std::shared_ptr<container_type> > m_steps;
        std::vector< scalar_type > m_times;
        std::vector< bool > m_out_of_bounds;
        bool m_verbose;
    };


    // tensor produce (structured) bspline field, for BARG
    // copied and modified from structured_mesh_time_dependent_field
    template<typename DataSet, typename T, unsigned int N, typename Container, typename Value = typename Container::value_type>
    class tp_bspline_time_dependent_field {
    public:
        typedef DataSet                                  dataset_type;  //this is boundaryAwareRectilinearGrid
        typedef vtk_utils::interpolator<DataSet, T, N>   interpolator_type;
        typedef typename interpolator_type::value_type   scalar_type;
        typedef typename interpolator_type::vector_type  point_type;
        typedef typename interpolator_type::vector_type  vector_type;
        typedef typename interpolator_type::matrix_type  matrix_type;
        typedef Container container_type;
        typedef Value value_type; //Eigen::Vector3d

    private:
        static constexpr scalar_type invalid_time = std::numeric_limits<scalar_type>::min();

        bool is_out_of_bounds(int cellid) const {
            if (m_out_of_bounds.empty()) return false;
            else return m_out_of_bounds[cellid];
        }

        value_type interpolate(int time_step,
            const std::vector<scalar_type>& weights,
            const std::vector<int>& ptids) const
        {
            std::shared_ptr<const container_type> values = m_steps[time_step];

            value_type v = value_type(weights.front()*(*values)[ptids.front()]);
            if (m_verbose) {
                std::cout << "ptid[0]=" << ptids[0] << ", weight[0]=" << weights.front()
                    << ", value[0]=" << to_str((*values)[ptids.front()]) << '\n';
            }
            for (int i = 1; i < weights.size(); i++) {
                v += value_type(weights[i] * (*values)[ptids[i]]);
                if (m_verbose) {
                    std::cout << "ptid[" << i << "]=" << ptids[i] << ", weight["
                        << i << "]=" << weights[i] << ", value[" << i << "]="
                        << to_str((*values)[ptids[i]]) << '\n';
                }
            }
            return v;
        }

    public:
        tp_bspline_time_dependent_field(bool verbose = false)
            : m_times(), m_steps(), m_verbose(verbose), m_out_of_bounds() {}

        tp_bspline_time_dependent_field(dataset_type* dataset,
            const std::vector< std::shared_ptr<container_type> >& steps,
            const std::vector< scalar_type >& times,
            const std::vector< bool >& out_of_bounds,
            bool verbose = false,
            bool use_bdry_aware = true)
            : m_interpolator(new interpolator_type(dataset)), m_steps(steps), m_times(times),
            m_verbose(verbose), m_out_of_bounds(out_of_bounds), m_dataset(dataset), m_bdry_aware(use_bdry_aware) {
            // check that the time steps are sorted and sort them if not
            if (!std::is_sorted(m_times.begin(), m_times.end())) {
                if (m_verbose) {
                    std::cout << "reordering the time steps" << std::endl;
                }

                std::vector<size_t> indices(m_times.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
                    return m_times[i1] < m_times[i2];
                });
                std::vector<scalar_type> new_times;
                std::vector< std::shared_ptr<container_type> > new_steps;
                for (size_t i = 0; i < indices.size(); ++i) {
                    new_times[i] = m_times[indices[i]];
                    new_steps[i] = m_steps[indices[i]];
                }
                m_times.swap(new_times);
                m_steps.swap(new_steps);
            }
        }

        std::pair<scalar_type, scalar_type> time_range() const {
            return std::pair<scalar_type, scalar_type>(m_times.front(), m_times.back());
        }

        value_type operator()(const point_type& p, scalar_type& t) const {
            // This is where the magic happens
            // this->m_dataset contains the boundaryAwareRectilinearGrid, and should be used for interpolation

            int idx = find_time(t, m_times);
            if (idx == -1) {
                std::ostringstream os;
                os << "tp_bspline_time_dependent_field::interpolate: invalid time coordinate: " << t << ", p=" << to_str(p);
                throw std::runtime_error(os.str());
            }
            else {
                // BARG version of interpolation. Note right now it only works for 1 time step.
                double query_param_tmp[3] = { p[0], p[1], p[2] };
                double return_val_tmp[3] = { 0 };
                bool is_starting_pt = t - m_times.front() < 1e-6;
                bool project_onto_closest_bdry_surface = m_bdry_aware && !is_starting_pt;
                vtkIdType ci = m_dataset->BsplineInterpolate(query_param_tmp, return_val_tmp, m_bdry_aware, project_onto_closest_bdry_surface);
                if (ci < 0) {
                    throw nvis::invalid_position_exception("BARG: out of bound.");
                }


                value_type ret_val(return_val_tmp[0], return_val_tmp[1], return_val_tmp[2]);
                return ret_val;

                /////// Below is the old version ////////

#if 0
                std::vector<scalar_type> weights(8);
                std::vector<int> ptids(8);
                int cellid;

                //interpolation_info gets the weights and the ptids for interpolation use
                if (!m_interpolator->interpolation_info(p, weights, ptids, cellid) ||
                    is_out_of_bounds(cellid)) {
                    std::ostringstream os;
                    os << "tp_bspline_time_dependent_field::interpolate: invalid location: " << to_str(p) << ", t=" << t;
                    if (m_verbose) {
                        std::cerr << os.str() << std::flush;
                    }
                    throw std::runtime_error(os.str());
                }

                value_type vlo = interpolate(idx, weights, ptids);      // spatial interpolation
                value_type vhi = interpolate(idx + 1, weights, ptids);    // spatial interpolation
                scalar_type u = (t - m_times[idx]) / (m_times[idx + 1] - m_times[idx]); // linear time interpolation

                if (m_verbose) {
                    std::ostringstream os;
                    os << "interpolation at " << to_str(p) << " at time step #" << idx << " returned " << to_str(vlo) << '\n';
                    os << "interpolation at " << to_str(p) << " at time step #" << idx + 1 << " returned " << to_str(vhi) << '\n';
                    std::cout << os.str() << std::flush;
                }
                vlo *= (1 - u);
                vhi *= u;
                return value_type(vlo + vhi);
#endif
            }
        }

        std::shared_ptr<interpolator_type> get_interpolator() {
            return m_interpolator;
        }

        vtkSmartPointer<dataset_type> m_dataset; //Added for interpolation calls
        std::shared_ptr< interpolator_type > m_interpolator;
        std::vector< std::shared_ptr<container_type> > m_steps;
        std::vector< scalar_type > m_times;
        std::vector< bool > m_out_of_bounds;
        bool m_verbose;
        bool m_bdry_aware;
    };


} // xavier

#endif
