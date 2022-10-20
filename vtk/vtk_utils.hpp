#ifndef __VTK_UTILS_HPP__
#define __VTK_UTILS_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace nvis {
    typedef bounding_box<vec1> bbox1;
}

#include <vtk/vtk_macros.hpp>
#include <vtk/vtk_data_helper.hpp>
#include <vtk/vtk_camera_helper.hpp>
#include <vtk/vtk_colorbar_helper.hpp>
#include <vtk/vtk_image_helper.hpp>
#include <vtk/vtk_io_helper.hpp>

#endif
