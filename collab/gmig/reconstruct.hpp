#ifndef _XAVIER_COLLAB_GMIG_SMOOTH_RECONSTRUCTION_HPP_
#define _XAVIER_COLLAB_GMIG_SMOOTH_RECONSTRUCTION_HPP_

// stl
#include <iostream>
#include <memory>
#include <vector>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>
// xavier
#include <data/raster.hpp>
#include "data_IO.hpp"
#include "utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xavier { namespace gmig {

/** \class travel_time_interpolator
  * \brief A lightweight wrapper around a RBF interpolator to solve
  * the travel time reconstruction problem.
  *
  * This class decomposes the travel time reconstruction problem into 
  * linear isotropic and nonlinear anisotropic part, solving the latter
  * with RBF interpolation. The former is determined by the average travel
  * velocity observed at all the receiving stations.
  *
  * \tparam BaseInt_  Type of the embedded RBF-based interpolating function
  */
template<typename BaseInt_>
class travel_time_interpolator {
public:
    typedef BaseInt_                             base_intp;
    typedef typename BaseInt_::function_type     rbf_type;
    typedef typename BaseInt_::point_type        point_type;
    typedef typename BaseInt_::value_type        data_type;
    typedef typename BaseInt_::scalar_type       scalar_type;
    typedef typename BaseInt_::derivative_type   derivative_type;
    typedef travel_time_data<scalar_type>        tt_data_type;
    
    travel_time_interpolator(const std::string& filename,
                             const rbf_type& phi = rbf_type(),
                             bool verbose = false);
    
    travel_time_interpolator(const tt_data_type& tt_data,
                             const rbf_type& phi = rbf_type(),
                             bool verbose=false);
                             
    data_type operator()(const point_type& x) const;
    derivative_type derivative(const point_type& x) const;
    
    const tt_data_type& data() const { return _tt_data; }

private:
    // helper functions
    void create_intp();
    
    std::shared_ptr<base_intp>  _intp;
    tt_data_type                _tt_data;
    rbf_type                    _phi;
    bool                        _verbose;
};

/** \fn check_solution
  * \brief Check interpolation quality
  *
  * Compute the mean and max relative error between a set of data points
  * and the values computed by the provided function at the corresponding
  * locations.
  *
  * \tparam Int_   Type of the interpolating function
  * \tparam Val_   Scalar value type for the reconstruction
  * \param fun     Interpolating function
  * \param points  Receivers' locations
  * \param times   Travel times measured at receivers
  * \return        Mean and maximum relative error (stored in a pair)
  */
template<typename Int_, typename Val_ = typename Int_::data_type>
std::pair<double, double>
check_solution(std::shared_ptr<const Int_> fun,
               const std::vector<nvis::fixed_vector<Val_, 2> >& points,
               const std::vector<nvis::fixed_vector<Val_, 1> >& times);

/** \fn travel_time
  * \brief Computes travel time at the vertices of a regular grid
  *
  * This function uses the provided function to determine the travel
  * time at the vertices of a regular grid.
  *
  * \tparam Int_   Type of the interpolating function
  * \tparam Val_   Scalar value type for the reconstruction
  * \param raster  Raster grid defining the vertices to be sampled and used
  *                to return the computed values
  * \param fun     Interpolating function
  * \param verbose Verbose flag (optional)
  */
template<typename Int_, typename Val_ = typename Int_::data_type>
void travel_time(raster2d<Val_>& raster, const Int_& fun,
                 bool verbose = false);

/** \fn travel_time
  * \brief Computes travel time and corresponding spatial gradient
  *        at the vertices of a regular grid
  *
  * This function uses the provided function to determine the travel
  * time at the vertices of a regular grid. In addition it computes the
  * spatial gradient of the travel time using the derivative of the function.
  *
  * \tparam Int_    Type of the interpolating function
  * \tparam Val_    Scalar value type for the reconstruction
  * \param  raster  Raster grid defining the vertices to be sampled and used
  *                 to return the computed values and derivatives
  * \param  fun     Interpolating function
  * \param  verbose Verbose flag (optional)
  */
template<typename Int_, typename Val_ = typename Int_::data_type>
void travel_time_gradient(raster2d<nvis::fixed_vector<Val_, 3> >& raster, 
                          const Int_& fun, bool verbose = false);
} // gmig
} // xavier


#include "detail/reconstruct.hpp"

#endif