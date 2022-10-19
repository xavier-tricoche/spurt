#ifndef __VTK_UTILS_HPP__
#define __VTK_UTILS_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace nvis {
    typedef bounding_box<vec1> bbox1;
}

#include <VTK/vtk_macros.hpp>
#include <VTK/vtk_data_helper.hpp>
#include <VTK/vtk_camera_helper.hpp>
#include <VTK/vtk_colorbar_helper.hpp>
#include <VTK/vtk_image_helper.hpp>
#include <VTK/vtk_io_helper.hpp>

#endif
