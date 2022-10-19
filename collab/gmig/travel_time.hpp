#ifndef __XAVIER_COLLAB_GMIG_TRAVEL_TIME_HPP__
#define __XAVIER_COLLAB_GMIG_TRAVEL_TIME_HPP__

#include <string>
#include <vector>

#include <math/fixed_vector.hpp>
#include "typedefs.hpp"
#include "utils.hpp"

namespace xavier { namespace gmig { namespace traveltime {
    
/** \struct travel_time_data
  * \brief Information needed to define the RBF based smooth reconstruction
  *        of travel time data points.
  *
  * This struct stores all the information needed to describe the input travel
  * time data and permit its smooth reconstruction using a combination of 
  * constant isotropic velocity and radial basis function interpolation or
  * approximation. The RBF kernel is identified by its name and may be 
  * accompanied by its radius if it is a compactly supported kernel.  
  *
  * \tparam Scalar_     Scalar type used to represent travel times and 
  *                     distances
  *
  * \var    source      Location of the source
  * \var    receivers   Locations of all the receivers
  * \var    distances   Distance of each receiver to source
  * \var    times       Travel times associated with the receivers
  * \var    velocity    Isotropic velocity component
  * \var    kernel_name Name of kernel used in RBF reconstruction
  * \var    radius      Radius of compactly supported RBF kernel
  * \var    weights     Weights of RBF reconstruction
  */
template<typename Scalar_ = double>
struct travel_time_data {
    typedef nvis::vec2 point;
    typedef Scalar_    scalar;
    
    travel_time_data() 
        : source(invalid_double), receivers(), times(), velocity(0), 
          kernel_name("none"), radius(0), weights() {}
          
    template<typename T>
    travel_time_data(const point& src, const std::vector<point>& recs,
                     const std::vector<T>& tt, 
                     const std::vector<T>& dist = std::vector<T>())
        : source(src), receivers(recs), times(tt), distances(dist) {
        if (!sanity_check()) {
            throw std::runtime_error("travel_time: sanity check failed!");
        }
        if (distances.empty()) {
            distances.resize(receivers.size());
            for (size_t i=0 ; i<receivers.size() ; ++i) {
                distances[i] = distance_in_km(source, receivers[i]);
            }
        }
        compute_velocity();
    }
          
    bool sanity_check() const {
        // check availability of original input data
        if (receivers.empty() || times.empty()) return false;
        if (receivers.size() != times.size()) return false;
        if (source[0] == invalid_double) return false;
        // check consistency of reconstruction information
        if (!weights.empty() && kernel_name == "none") return false;
        return true;
    }
    
    double compute_velocity() {
        double _sum = 0;
        for (size_t i=0 ; i<receivers.size() ; ++i) {
            double dist = distance_in_km(receivers[i], source);
            _sum += dist / times[i];
        }
        velocity = _sum / (double)receivers.size();
        return velocity;
    }
    
    void clear() {
        source[0] = invalid_double;
        receivers.clear();
        times.clear();
        distances.clear();
        velocity = 0;
        kernel_name = "none";
        radius = 0;
        weights.clear();
    }
    
    size_t size() const {
        return receivers.size();
    }
    
    // original data
    point               source;
    std::vector<point>  receivers;
    std::vector<scalar> times;
    std::vector<scalar> distances;
    
    // reconstruction data
    double              velocity;
    std::string         kernel_name;
    double              radius;
    std::vector<scalar> weights;
};

template<typename Int_, typename Val_ = typename Int_::scalar_type>
struct interpolator {
    typedef Int_                          rbf_type;
    typedef Val_                          scalar;
    typedef travel_time_data<scalar>      tt_data_type;
    typedef typename tt_data_type::point  point;
    
    interpolator(const tt_data_type& ttdata, bool verbose=false,
                 const typename Int_::function_type& phi = 
                       typename Int_::function_type())
        : _source(ttdata.source), _velocity(ttdata.velocity) {
        std::vector<scalar> aniso_times(ttdata.times.begin(),
                                        ttdata.times.end());
        for (size_t i=0 ; i<aniso_times.size() ; ++i) {
            aniso_times[i] -= ttdata.distances[i]/_velocity;
        }
        _rbf = rbf_type(ttdata.receivers, aniso_times, phi, verbose);
    }
    
    double operator()(const point& x) const {
        return distance_in_km(x, _source)/_velocity + _rbf(x);
    }
    
    rbf_type  _rbf;
    point     _source;
    double    _velocity;
};

} // namespace traveltime
} // namespace gmig
} // namespace xavier

#endif