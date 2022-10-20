#ifndef __XAVIER_COLLAB_MAARTEN_CORE_MANTEL_IO_HPP__
#define __XAVIER_COLLAB_MAARTEN_CORE_MANTEL_IO_HPP__

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
// spurt
#include "typedefs.hpp"
#include "core_mantle.hpp"

namespace spurt { namespace gmig { namespace core_mantle {
    
/** \fn read_text
  * \brief Custom reader for core mantle text files (XYZ.txt)
  * 
  * Import the location of points corresponding to various
  * interfaces around the core-mantle boundary. The 3D position of 
  * each point is imported along with the associated layer ID and type. 
  *
  * \code
  * format:
  * layer_no; lat; lon; height; type; lat_sec; lon_sec.
  * \endcode
  * 
  * \tparam Scalar_   Scalar type used to represent position coordinates
  * \param  data      Imported interface dataset
  * \param  file_name File to read from
  * \param  verbose   Verbose flag (optional)
  * \return           Bouding box of receivers' locations
  */
template<typename Scalar_ = nvis::vec1, typename Index_ = int>
typename Vertex<Scalar_,Index_>::bbox_t 
read_text(std::vector<Vertex<Scalar_,Index_> >& data,
          const std::string& file_name,
          bool verbose = false);
} // namespace core_mantle
} // namespace gmig
} // namespace spurt

#include "detail/core_mantle_io.hpp"

#endif