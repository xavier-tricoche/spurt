#ifndef __XAVIER_COLLAB_MAARTEN_DATA_IO_HPP__
#define __XAVIER_COLLAB_MAARTEN_DATA_IO_HPP__

// STL
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
// xavier
#include "typedefs.hpp"

namespace xavier { namespace gmig { namespace traveltime {

/** \fn read_text
  * \brief Custom reader for Hui Huang's files
  * 
  * Import the location of the receivers and the associated value
  * of the travel time along with the location of the source, which
  * is returned as the first position with travel time 0. If available, 
  * the explicit distance between source and each receiver is imported
  * as well. The format is as follows:
  *
  * \code
  * format:
  * source_lon; source_lat; receiver_lon; receiver_lat; travel_time; \
  * <optional: source-receiver distance>; receiver index.
  * \endcode
  * 
  * \tparam Scalar_   Scalar type used to represent travel times and distances
  * \param  data      Imported travel time dataset
  * \param  file_name File to read from
  * \param  verbose   Verbose flag (optional)
  * \return           Bouding box of receivers' locations
  */
template<typename Scalar_ = nvis::vec1>
nvis::bbox2 read_text(travel_time_data<Scalar_>& data,
                      const std::string& file_name,
                      bool verbose = false);
                      
/** \fn read_nrrd
  * \brief Import travel time data and RBF reconstruction information
  *        stored in NRRD format
  * 
  * Import the location of the receivers and the associated travel time value.
  * If the parameters of the RBF reconstruction were previously computed and 
  * stored, they are imported as well. The expected format is the one produced 
  * by save_rbf().
  * 
  * \tparam Scalar_   Scalar type used to represent travel times and distances
  * \param  data      Imported travel time dataset
  * \param  file_name File to read from
  * \param  verbose   Verbose flag (optional)
  * \return           Bouding box of receivers' locations
  */
template<typename Scalar_ = nvis::vec1>
nvis::bbox2
read_nrrd(travel_time_data<Scalar_>& data,
          const std::string& filename,
          bool verbose = false);
          
/** \fn save_rbf
  * \brief Export RBF reconstruction data in NRRD format
  * 
  * Exports in NRRD format all the information necessary to smoothly
  * reconstruct the travel time information. This format 
  * is the one read by read_nrrd().
  * 
  * \tparam Scalar_   Scalar type used to represent travel times and distances
  * \param  data      Travel time dataset and reconstruction information
  * \param  file_name File to read
  * \param  verbose   Verbose flag (optional)
  */
template<typename Scalar_ = nvis::vec1>
void save_rbf(const travel_time_data<Scalar_>& data,
              const std::string& filename,
              bool verbose = false);
              
} // namespace traveltime
} // namespace gmig
} // namespace xavier
                      
#include "detail/data_IO.hpp"

#endif