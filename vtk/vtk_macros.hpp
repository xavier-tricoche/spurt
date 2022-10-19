#ifndef __VTK_MACROS_HPP__
#define __VTK_MACROS_HPP__

#include <vtkAlgorithm.h>
#include <vtkCharArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkGraph.h>
#include <vtkLongArray.h>
#include <vtkLongLongArray.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRectilinearGrid.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredPoints.h>
#include <vtkUniformGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkUnstructuredGrid.h>

#include <misc/meta_utils.hpp>

// Smart pointer associated with a VTK type
#define VTK_SMART(type) \
    vtkSmartPointer<type>

// Declare and initialize a VTK pointer
#define VTK_PTR(type, name)  \
    type* name = type::New();

// Declare and initialize a VTK smart pointer
#define VTK_CREATE(type, name)  \
    vtkSmartPointer<type> name = vtkSmartPointer<type>::New();
    
// Initialize an existing VTK smart pointer
#define VTK_INIT(type, name) \
    name = vtkSmartPointer<type>::New();

// Connect output ports and input connections
#define VTK_PLUG(receiver, provider)  \
    (receiver)->SetInputConnection((provider)->GetOutputPort());

// Connect output ports and input source
#define VTK_SOURCE_PLUG(receiver, provider) \
    (receiver)->SetSourceConnection((provider)->GetOutputPort());

// Same but for regular input / output pointers
// Fix for backward incompatibility between VTK6 and VTK5 required here
#if (VTK_MAJOR_VERSION >= 6)
#define VTK_CONNECT(receives, supplies) \
    (receives)->SetInputData(supplies);
#else
#define VTK_CONNECT(receives, supplies) \
    (receives)->SetInput(supplies);
#endif

// Special case of Image to Texture connections
#if (VTK_MAJOR_VERSION >= 6)
#define VTK_IMAGE_CONNECT(receives, image) \
    (receives)->SetInputData(image);
#else
#define VTK_IMAGE_CONNECT(receives, image) \
    (receives)->SetInputConnection((image)->GetProducerPort());
#endif

// Special case of Source connections
#if (VTK_MAJOR_VERSION >= 6)
#define VTK_SOURCE_CONNECT(receives, supplies) \
    (receives)->SetSourceData(supplies);
#else
#define VTK_SOURCE_CONNECT(receives, supplies) \
    (receives)->SetSource(supplies);
#endif
    
// Create an actor pointer for a given object.
// NOTE: A mapper of the correct type is automatically
//       deduced and created. It is accessible upon 
//       completion of this macro through the actor's
//       GetMapper() method.
#define VTK_MAKE_ACTOR(actor, object) \
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New(); \
    actor->SetMapper(vtk_utils::right_mapper(vtk_utils::get_ptr(object)));

namespace vtk_utils {

// X-Macro used to automate the creation of type traits for vtkDataObject subtypes
#define LIST_OF_MAPABLE_VTK_DATA_OBJECT_TYPES \
    X(vtkPolyData, VTK_POLY_DATA, vtkPolyDataMapper); \
    X(vtkStructuredPoints, VTK_STRUCTURED_POINTS, vtkDataSetMapper); \
    X(vtkStructuredGrid, VTK_STRUCTURED_GRID, vtkDataSetMapper); \
    X(vtkRectilinearGrid, VTK_RECTILINEAR_GRID, vtkDataSetMapper); \
    X(vtkUnstructuredGrid, VTK_UNSTRUCTURED_GRID, vtkDataSetMapper); \
    X(vtkDataSet, VTK_DATA_SET, vtkDataSetMapper); \
    X(vtkPointSet, VTK_POINT_SET, vtkDataSetMapper); \
    X(vtkUniformGrid, VTK_UNIFORM_GRID, vtkDataSetMapper); \
    
template<typename T>
struct vtk_mapable_data_object_traits {};

#define X(type, typeID, mapper) \
    template<> \
    struct vtk_mapable_data_object_traits<type> { \
        typedef type   data_object_type; \
        typedef mapper mapper_type; \
    };
    LIST_OF_MAPABLE_VTK_DATA_OBJECT_TYPES
#undef X

#define LIST_OF_VTK_DATA_TYPES \
    X(char, Char); \
    X(unsigned char, UnsignedChar); \
    X(short, Short); \
    X(unsigned short, UnsignedShort); \
    X(int, Int); \
    X(unsigned int, UnsignedInt); \
    X(long, Long); \
    X(unsigned long, UnsignedLong); \
    X(long long, LongLong); \
    X(unsigned long long, UnsignedLongLong); \
    X(float, Float); \
    X(double, Double);

template<typename T>
struct vtk_array_traits {};

#define X(type, typestr) \
    template<> \
    struct vtk_array_traits<type> { \
        typedef type                data_type; \
        typedef vtk##typestr##Array array_type; \
    };
    LIST_OF_VTK_DATA_TYPES
#undef X

template< typename DataObject, 
          typename Mapper = 
          typename vtk_mapable_data_object_traits<DataObject>::mapper_type >
inline Mapper* right_mapper(DataObject* data) {
    Mapper* mapper = Mapper::New();
#if VTK_MAJOR_VERSION >= 6
    mapper->SetInputData(data);
#else
    mapper->SetInput(data);
#endif
    return mapper;
}

template<typename T>
inline T* get_ptr(vtkSmartPointer<T> sptr) {
    return sptr.GetPointer();
}

template<typename T>
inline T* get_ptr(T* ptr) {
    return ptr;
}

inline void vtk_connect(vtkAlgorithm* receiver, vtkAlgorithm* provider) {
    receiver->SetInputConnection(provider->GetOutputPort());
}

inline void vtk_connect(vtkAlgorithm* receiver, vtkDataObject* provider) {
#if VTK_MAJOR_VERSION >= 6    
    receiver->SetInputDataObject(provider);
#else
    receiver->SetInputConnection(provider->GetProducerPort());
#endif
}

} // namespace vtk_utils

#endif
