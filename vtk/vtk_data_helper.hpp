#ifndef __VTK_DATA_HELPER_HPP__
#define __VTK_DATA_HELPER_HPP__

#include "vtk_macros.hpp"

#include "vtkActor.h"
#include "vtkAlgorithm.h"
#include "vtkBitArray.h"
#include "vtkBMPReader.h"
#include "vtkBMPWriter.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLinks.h"
#include "vtkCharArray.h"
#include "vtkClipPolyData.h"
#include "vtkColorTransferFunction.h"
#include "vtkCommand.h"
#include "vtkContourFilter.h"
#include "vtkCutter.h"
#include "vtkCylinderSource.h"
// #include "vtkDataArrayTemplate.h"
#include "vtkDataSetMapper.h"
#include "vtkDataSetReader.h"
#include "vtkDataSetWriter.h"
#include "vtkDelaunay2D.h"
#include "vtkDelaunay3D.h"
#include "vtkDoubleArray.h"
#include "vtkExtractEdges.h"
#include "vtkFloatArray.h"
#include "vtkGraph.h"
#include "vtkGlyph2D.h"
#include "vtkGlyph3D.h"
#include "vtkGlyphSource2D.h"
#include "vtkImageCast.h"
#include "vtkImageData.h"
// #include "vtkImageDataLIC2D.h"
#include "vtkImageShiftScale.h"
#include "vtkIntArray.h"
#include "vtkInteractorObserver.h"
#include "vtkJPEGReader.h"
#include "vtkJPEGWriter.h"
#include "vtkLightKit.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#if VTK_MAJOR_VERSION >= 6
#include "vtkNrrdReader.h"
#endif
#include "vtkObject.h"
#include "vtkPlane.h"
#include "vtkPlaneSource.h"
#include "vtkPNGReader.h"
#include "vtkPNGWriter.h"
#include "vtkPNMReader.h"
#include "vtkPNMWriter.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataReader.h"
#include "vtkProbeFilter.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkShortArray.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkSphereSource.h"
#include "vtkStaticCellLinks.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkTIFFReader.h"
#include "vtkTIFFWriter.h"
#include "vtkTransform.h"
#include "vtkTubeFilter.h"
#include "vtkUniformGrid.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkWindowToImageFilter.h"

#if __VTK_HAS_DAG__
#include "vtkDirectedGraph.h"
#include "vtkDirectedAcyclicGraph.h"
#include "vtkGraphMapper.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkUndirectedGraph.h"
#endif

#include <cassert>
#include <list>
#include <set>
#include <string>
#include <stdexcept>
#include <vector>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
//#include <type_traits>
#include <format/filename.hpp>
#include <image/nrrd_wrapper.hpp>
#include <math/vector_manip.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

namespace nvis {
    typedef fixed_vector<double, 1> vec1;
    typedef bounding_box<vec1>     bbox1;
}


namespace vtk_utils {

template<typename PointOut, typename PointIn>
inline PointOut make3d(const PointIn& _in)
{
    const size_t N = _in.size();
    assert(N>0);
    switch (N) {
        case 1:
            return PointOut(_in[0], 0, 0);
        case 2:
            return PointOut(_in[0], _in[1], 0);
        default:
            return PointOut(_in[0], _in[1], _in[2]);
    }
}

template<typename DataSetPtr>
inline DataSetPtr add_vertices(DataSetPtr);

template< typename ForwardContainer>
inline vtkPoints* make_vtkpoints(const ForwardContainer& pos)
{
    typedef typename ForwardContainer::value_type       point_type;
    typedef typename ForwardContainer::const_iterator   iterator_type;
    typedef typename point_type::value_type                    value_type;
    typedef typename vtk_array_traits<value_type>::array_type  array_type;
    typedef typename nvis::fixed_vector<value_type, 3>         vec_type;

    VTK_CREATE(array_type, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(pos.size());
    vtkIdType _count=0;
    for (iterator_type it=pos.begin() ; it!=pos.end() ; ++it, ++_count) {
        vec_type x = make3d<vec_type, point_type>(*it);
        coords->SetTypedTuple(_count, &x[0]);
    }
    VTK_PTR(vtkPoints, points);
    points->SetData(coords);
    return points;
}

// Return a polydata that contains only point information
template< typename ForwardContainer>
inline vtkPolyData* make_points(const ForwardContainer& pos)
{
    vtkPoints* points = make_vtkpoints(pos);
    VTK_PTR(vtkPolyData, polydata);
    polydata->SetPoints(points);
    return polydata;
}

// Return a polydata that contains only point information
// Here the points correspond to a subset of the input array of points
// as defined by the mask "selected"
template< typename RandomAccessContainer,
          typename ForwardContainer >
inline vtkPolyData* make_points(const RandomAccessContainer& pos,
                                const ForwardContainer& selected)
{
    typedef typename RandomAccessContainer::value_type         point_type;
    typedef typename ForwardContainer::const_iterator          iterator_type;
    typedef typename point_type::value_type                    value_type;
    typedef typename vtk_array_traits<value_type>::array_type  array_type;
    typedef typename nvis::fixed_vector<value_type, 3>         vec_type;

    VTK_CREATE(array_type, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(selected.size());
    vtkIdType _count=0;
    for (iterator_type it=selected.begin() ; it!=selected.end() ; ++it, ++_count) {
        vec_type x = make3d<vec_type, point_type>(pos[*it]);
        coords->SetTypedTuple(_count, &x[0]);
    }
    VTK_CREATE(vtkPoints, points);
    points->SetData(coords);
    VTK_PTR(vtkPolyData, polydata);
    polydata->SetPoints(points);
    return polydata;
}

// Return a polydata that contains polylines
template< typename ForwardContainer>
inline vtkPolyData*
make_polylines(const ForwardContainer& lines,
               std::vector<nvis::ivec2>& removed, double mind=1.0e-6,
               bool verbose=false)
{
    typedef typename ForwardContainer::const_iterator          iterator_type;
    typedef typename ForwardContainer::value_type              line_type;
    typedef typename line_type::const_iterator                 line_iterator_type;
    typedef typename line_type::value_type                     point_type;
    typedef typename point_type::value_type                    value_type;
    typedef typename vtk_array_traits<value_type>::array_type  array_type;
    typedef typename nvis::fixed_vector<value_type, 3>         vec_type;

    removed.clear();

    VTK_CREATE(array_type, coords);
    coords->SetNumberOfComponents(3);
    VTK_CREATE(vtkCellArray, cells);

    size_t lid=0;
    for (iterator_type it=lines.begin(); it!=lines.end(); ++it, ++lid) {
        const line_type& line = *it;
        if (line.size() < 2) {
            std::cerr << "Invalid number of points in line #" << lid << ": " << line.size() << '\n';
            continue;
        }

        // run sanity check on line
        std::vector<point_type> valid_pts;
        int i=0;
        for (line_iterator_type lt=line.begin() ; lt!=line.end() ; ++lt, ++i) {
            if (valid_pts.empty()) {
                valid_pts.push_back(*lt);
            }
            else {
                const point_type& p=*lt;
                const point_type& q=valid_pts.back();
                if (!mind || // no lower bound on segment length
                    (i==line.size()-1) || // always include last point
                    (mind>0 && spurt::vector::distance(p, q)>mind)) { // else we must pass the filtering criterion
                    valid_pts.push_back(p);
                }
                else {
                    removed.push_back(nvis::ivec2(lid, i));
                    if (verbose) {
                        std::cout << "point #" << i << " rejected, distance="
                            << spurt::vector::distance(p, q) << '\n';
                        std::cout << "p=" << p << ", q=" << q << '\n';
                    }
                }
            }
        }

        cells->InsertNextCell(valid_pts.size());
        for (auto it=valid_pts.begin(); it!=valid_pts.end(); ++it) {
            const point_type& p = *it;
            vec_type x = make3d<vec_type, point_type>(p);
            cells->InsertCellPoint(coords->InsertNextTypedTuple(&x[0]));
        }
    }
    VTK_CREATE(vtkPoints, points);
    points->SetData(coords);
    VTK_PTR(vtkPolyData, polydata);
    polydata->SetPoints(points);
    polydata->SetLines(cells);
    return polydata;
}

// Return a polydata that contains polylines
template< typename ForwardContainer>
inline vtkPolyData*
make_polylines(const ForwardContainer& lines, double mind=1.0e-6,
               bool verbose=false)
{
    std::vector<nvis::ivec2> dummy;
    return make_polylines(lines, dummy, mind, verbose);
}

// Return a polydata that contains polygons
template< typename ForwardContainer>
inline vtkPolyData* add_edges_from_numbers(vtkPolyData* inout,
        const ForwardContainer& edges)
{
    VTK_CREATE(vtkCellArray, cells);
    for (auto it=edges.begin(); it!=edges.end(); ++it) {
        cells->InsertNextCell(2);
        cells->InsertCellPoint(*it); ++it;
        cells->InsertCellPoint(*it);
    }
    inout->SetLines(cells);
    return inout;
}

template< typename DataSetPtr >
inline vtkPolyData* make_spheres(DataSetPtr in, double radius,
                                int theta_res=12, int phi_res=12,
                                bool scale=false)
{
    VTK_PTR(vtkSphereSource, source);
    source->SetRadius(radius);
    source->SetPhiResolution(phi_res);
    source->SetThetaResolution(theta_res);
    VTK_PTR(vtkGlyph3D, glyphs);
    VTK_CONNECT(glyphs, in);
    glyphs->SetSourceConnection(source->GetOutputPort());
    glyphs->SetScaling(scale);
    glyphs->SetScaleModeToScaleByScalar();
    glyphs->Update();
    vtkPolyData *spheres = glyphs->GetOutput();
    spheres->Register(0);

    glyphs->Delete();
    source->Delete();

    return spheres;
}

// add scalar attributes to point/cell data
template< typename ForwardContainer, typename DataSetPtr >
inline DataSetPtr add_scalars(DataSetPtr inout, const ForwardContainer& scalars,
                              bool point_data=true, const std::string& name="anonymous_scalars",
                              bool active=true)
{
    typedef typename ForwardContainer::value_type             value_type;
    typedef typename ForwardContainer::const_iterator         iterator_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;

    VTK_CREATE(array_type, values);
    values->SetNumberOfComponents(1);
    values->SetNumberOfTuples(scalars.size());
    vtkIdType _count=0;
    for (iterator_type it=scalars.begin() ; it!=scalars.end() ; ++it, ++_count) {
        values->SetTypedTuple(_count, &*it);
    }
    values->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetScalars(values);
        else inout->GetPointData()->AddArray(values);
    }
    else {
        if (active) inout->GetCellData()->SetScalars(values);
        else inout->GetCellData()->AddArray(values);
    }
    return inout;
}

template< typename T, typename DataSetPtr>
inline DataSetPtr add_scalars_from_carray(DataSetPtr inout, const T* scalars, size_t nvals,
                                          bool point_data=true, const std::string& name="anonymous_scalars",
                                          bool active=true) {
    typedef T value_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;
    VTK_CREATE(array_type, values);
    values->SetNumberOfComponents(1);
    values->SetNumberOfTuples(nvals);
    for (size_t i=0; i<nvals; ++i) {
        values->SetTypedTuple(i, &scalars[i]);
    }
    values->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetScalars(values);
        else inout->GetPointData()->AddArray(values);
    }
    else {
        if (active) inout->GetCellData()->SetScalars(values);
        else inout->GetCellData()->AddArray(values);
    }
    return inout;
}

// add color attributes to point/cell data. "colors" is a container of RGB
// 3-vectors that can be cast to unsigned char
template< typename ForwardContainer, typename DataSetPtr >
inline DataSetPtr add_colors(DataSetPtr inout, const ForwardContainer& colors,
                             bool point_data=true, const std::string& name="anonymous_colors",
                             bool active=true)
{
    typedef typename ForwardContainer::value_type color_type;
    typedef typename ForwardContainer::const_iterator iterator_type;
    typedef typename vtk_array_traits<unsigned char>::array_type array_type;

    VTK_CREATE(array_type, values);
    values->SetNumberOfComponents(3);
    values->SetNumberOfTuples(colors.size());
    unsigned char _col[3];
    vtkIdType _count=0;
    for (iterator_type it=colors.begin() ; it!=colors.end() ; ++it, ++_count) {
        _col[0] = static_cast<unsigned char>( (*it)[0] );
        _col[1] = static_cast<unsigned char>( (*it)[1] );
        _col[2] = static_cast<unsigned char>( (*it)[2] );
        values->SetTypedTuple(_count, _col);
    }
    values->SetName(name.c_str());
    // special case of explicit colors passed on as scalars
    if (point_data) {
        if (active) inout->GetPointData()->SetScalars(values);
        else inout->GetPointData()->AddArray(values);
    }
    else {
        if (active) inout->GetCellData()->SetScalars(values);
        else inout->GetCellData()->AddArray(values);
    }
    return inout;
}

// Add vector attributes to point/cell data. "vectors" is a container
// of 1D arrays
template< typename ForwardContainer,
          typename DataSetPtr >
inline DataSetPtr add_vectors(DataSetPtr inout, const ForwardContainer& vectors,
                              bool point_data=true, const std::string& name="anonymous_vectors",
                              bool active=true)
{
    typedef typename ForwardContainer::value_type             vector_type;
    typedef typename ForwardContainer::const_iterator         iterator_type;
    typedef typename vector_type::value_type                  value_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;

    size_t N=0;
    if (!vectors.empty()) N = vectors.front().size();

    assert(N==2 || N==3);

    VTK_CREATE(array_type, _vectors);
    _vectors->SetNumberOfComponents(3);
    _vectors->SetNumberOfTuples(vectors.size());
    vtkIdType _count = 0;
    value_type _vec[3] = {0, 0, 0};
    if (N==2) {
        for (iterator_type it=vectors.begin() ; it!=vectors.end() ; ++it, ++_count) {
            _vec[0] = (*it)[0];
            _vec[1] = (*it)[1];
            _vectors->SetTypedTuple(_count, _vec);
        }
    } else {
        for (iterator_type it=vectors.begin() ; it!=vectors.end() ; ++it, ++_count) {
            _vec[0] = (*it)[0];
            _vec[1] = (*it)[1];
            _vec[2] = (*it)[2];
            _vectors->SetTypedTuple(_count, _vec);
        }
    }
    _vectors->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetVectors(_vectors);
        else inout->GetPointData()->AddArray(_vectors);
    }
    else {
        if (active) inout->GetCellData()->SetVectors(_vectors);
        else inout->GetCellData()->AddArray(_vectors);
    }
    return inout;
}

// Add vector attributes to point/cell data. "vectors" is a container of
// vector components, each of size N
template< typename ForwardContainer, typename DataSetPtr, int N >
inline DataSetPtr add_vectors_from_numbers(DataSetPtr inout,
                                           const ForwardContainer& vectors,
                                           bool point_data=true,
                                           const std::string& name="anonymous_vectors",
                                           bool active=true)
{
    static_assert(N==2 || N==3, "invalid vector dimensions in add_vectors_from_numbers");

    typedef typename ForwardContainer::const_iterator         iterator_type;
    typedef typename ForwardContainer::value_type             value_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;

    VTK_CREATE(array_type, _vectors);
    _vectors->SetNumberOfComponents(3);
    _vectors->SetNumberOfTuples(vectors.size()/N);
    vtkIdType _count = 0;
    value_type _vec[3] = {0, 0, 0};
    if (N==2) {
        for (iterator_type it=vectors.begin() ; it!=vectors.end() ; ++_count) {
            _vec[0] = (*it++);
            _vec[1] = (*it++);
            _vectors->SetTypedTuple(_count, _vec);
        }
    } else {
        for (iterator_type it=vectors.begin() ; it!=vectors.end() ; ++_count) {
            _vec[0] = (*it++);
            _vec[1] = (*it++);
            _vec[2] = (*it++);
            _vectors->SetTypedTuple(_count, _vec);
        }
    }
    _vectors->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetVectors(_vectors);
        else inout->GetPointData()->AddArray(_vectors);
    }
    else {
        if (active) inout->GetCellData()->SetVectors(_vectors);
        else inout->GetCellData()->AddArray(_vectors);
    }
    return inout;
}

template< typename T, typename DataSetPtr, int N>
inline DataSetPtr add_vectors_from_carray(DataSetPtr inout, const T* vectors, size_t nvals,
                                          bool point_data=true, const std::string& name="anonymous_vectors",
                                          bool active=true) {
    typedef T value_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;
    VTK_CREATE(array_type, _vectors);
    _vectors->SetNumberOfComponents(3);
    _vectors->SetNumberOfTuples(nvals);
    value_type v[3] = {0, 0, 0};
    if (N==2) {
        for (size_t i=0; i<nvals; ++i) {
            v[0] = vectors[2*i];
            v[1] = vectors[2*i+1];
            _vectors->SetTypedTuple(i, v);
        }
    }
    else {
        for (size_t i=0; i<nvals; ++i) {
            v[0] = vectors[3*i];
            v[1] = vectors[3*i+1];
            v[2] = vectors[3*i+2];
            _vectors->SetTypedTuple(i, v);
        }
    }
    _vectors->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetVectors(_vectors);
        else inout->GetPointData()->AddArray(_vectors);
    }
    else {
        if (active) inout->GetCellData()->SetVectors(_vectors);
        else inout->GetCellData()->AddArray(_vectors);
    }
    return inout;
}

// Add tensor attributes to point/cell data. "tensors" is a container
// of 1D arrays.
template< typename ForwardContainer, typename DataSetPtr >
inline DataSetPtr add_tensors(DataSetPtr inout, const ForwardContainer& tensors,
                              bool point_data=true, const std::string& name="anonymous_tensors",
                              bool active=true)
{
    typedef typename ForwardContainer::value_type             vector_type;
    typedef typename ForwardContainer::const_iterator         iterator_type;
    typedef typename vector_type::value_type                  value_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;
    static const size_t N = vector_type::size();

    assert(N==3 || N==4 || N==6 || N==9);

    VTK_CREATE(array_type, _tensors);
    _tensors->SetNumberOfComponents(9);
    _tensors->SetNumberOfTuples(tensors.size()/N);
    vtkIdType _count = 0;
    value_type _ten[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    if (N==3) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            _ten[0] = (*it++);
            _ten[1] = (*it++);
            _ten[3] = _ten[1];
            _ten[4] = (*it++);
            _tensors->SetTypedTuple(_count, _ten);
        }
    } else if (N==4) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            _ten[0] = (*it++);
            _ten[1] = (*it++);
            _ten[3] = (*it++);
            _ten[4] = (*it++);
            _tensors->SetTypedTuple(_count, _ten);
        }
    } else if (N==6) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            _ten[0] = (*it++);
            _ten[1] = (*it++);
            _ten[2] = (*it++);
            _ten[3] = _ten[1];
            _ten[4] = (*it++);
            _ten[5] = (*it++);
            _ten[6] = _ten[2];
            _ten[7] = _ten[5];
            _ten[8] = (*it++);
            _tensors->SetTypedTuple(_count, _ten);
        }
    } else if (N==9) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            for (size_t i=0 ; i<9 ; ++i) {
                _ten[i] = (*it++);
            }
            _tensors->SetTypedTuple(_count, _ten);
        }
    }
    _tensors->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetTensors(_tensors);
        else inout->GetPointData()->AddArray(_tensors);
    }
    else {
        if (active) inout->GetCellData()->SetTensors(_tensors);
        else inout->GetCellData()->AddArray(_tensors);
    }
    return inout;
}

// Add tensor attributes to point/cell data. "tensors" is a container
// of tensor components, each of size N.
template< typename ForwardContainer, typename DataSetPtr, size_t N >
inline DataSetPtr add_tensors_from_numbers(DataSetPtr inout,
                                           const ForwardContainer& tensors,
                                           bool point_data=true,
                                           const std::string& name="anonymous_tensors",
                                           bool active=true)
{
    static_assert(N==3 || N==4 || N==6 || N==9, "invalid tensor dimensions in add_tensors_from_numbers");

    typedef typename ForwardContainer::value_type             value_type;
    typedef typename ForwardContainer::const_iterator         iterator_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;

    VTK_CREATE(array_type, _tensors);
    _tensors->SetNumberOfComponents(9);
    _tensors->SetNumberOfTuples(tensors.size()/N);
    vtkIdType _count = 0;
    value_type _ten[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    if (N==3) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            _ten[0] = (*it++);
            _ten[1] = (*it++);
            _ten[3] = _ten[1];
            _ten[4] = (*it++);
            _tensors->SetTypedTuple(_count, _ten);
        }
    } else if (N==4) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            _ten[0] = (*it++);
            _ten[1] = (*it++);
            _ten[3] = (*it++);
            _ten[4] = (*it++);
            _tensors->SetTypedTuple(_count, _ten);
        }
    } else if (N==6) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            _ten[0] = (*it++);
            _ten[1] = (*it++);
            _ten[2] = (*it++);
            _ten[3] = _ten[1];
            _ten[4] = (*it++);
            _ten[5] = (*it++);
            _ten[6] = _ten[2];
            _ten[7] = _ten[5];
            _ten[8] = (*it++);
            _tensors->SetTypedTuple(_count, _ten);
        }
    } else if (N==9) {
        for (iterator_type it=tensors.begin() ; it!=tensors.end() ; ++_count) {
            for (size_t i=0 ; i<9 ; ++i) {
                _ten[i] = (*it++);
            }
            _tensors->SetTypedTuple(_count, _ten);
        }
    }
    _tensors->SetName(name.c_str());
    if (point_data) {
        if (active) inout->GetPointData()->SetTensors(_tensors);
        else inout->GetPointData()->AddArray(_tensors);
    }
    else {
        if (active) inout->GetCellData()->SetTensors(_tensors);
        else inout->GetCellData()->AddArray(_tensors);
    }
    return inout;
}

// Add texture coordinates to point data. "tcoords" is a container
// of 1D arrays
template< typename ForwardContainer,
          typename DataSetPtr >
inline DataSetPtr add_tcoords(DataSetPtr inout, const ForwardContainer& tcoords)
{
    typedef typename ForwardContainer::value_type             coord_type;
    typedef typename ForwardContainer::const_iterator         iterator_type;
    typedef typename coord_type::value_type                   value_type;
    typedef typename vtk_array_traits<value_type>::array_type array_type;

    VTK_CREATE(array_type, _tcoords);
    _tcoords->SetNumberOfComponents(2);
    _tcoords->SetNumberOfTuples(tcoords.size());
    vtkIdType _count = 0;
    value_type _coord[2] = {0, 0};
    for (iterator_type it=tcoords.begin() ; it!=tcoords.end() ; ++it, ++_count) {
        _coord[0] = (*it)[0];
        _coord[1] = (*it)[1];
        _tcoords->SetTypedTuple(_count, _coord);
    }
    inout->GetPointData()->SetTCoords(_tcoords);
    return inout;
}

template<typename DataSetPtr>
inline DataSetPtr add_vertices(DataSetPtr inout)
{
    VTK_CREATE(vtkCellArray, vertices);
    size_t npts = inout->GetNumberOfPoints();
    vertices->InitTraversal();
    for (size_t n=0 ; n<npts ; ++n) {
        vertices->InsertNextCell(1);
        vertices->InsertCellPoint(n);
    }
    inout->SetVerts(vertices);
    return inout;
}

template< typename DataSetPtr,
          typename ForwardContainer >
inline DataSetPtr add_lines(DataSetPtr inout, const ForwardContainer& lines)
{
    typedef typename ForwardContainer::const_iterator iterator_type;

    VTK_CREATE(vtkCellArray, cells);
    for (iterator_type it=lines.begin(); it!=lines.end(); ) {
        cells->InsertNextCell(2);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    inout->SetLines(cells);
    return inout;
}

template< typename DataSetPtr,
          typename ForwardContainer >
inline DataSetPtr add_polylines(DataSetPtr inout, const ForwardContainer& lines)
{
    typedef typename ForwardContainer::value_type line_type;
    typedef typename line_type::value_type id_type;

    VTK_CREATE(vtkCellArray, cells);
    for (auto it=lines.begin(); it!=lines.end(); ++it) {
        cells->InsertNextCell(it->size());
        std::for_each(it->begin(), it->end(), [&](const id_type& id)
            {
                cells->InsertCellPoint(id);
            });
    }
    inout->SetLines(cells);
    return inout;
}

inline vtkSmartPointer<vtkPolyData>
add_mesh2d(vtkSmartPointer<vtkPolyData> inout)
{

#if 0
    VTK_CREATE(vtkPolyData, copy);
    copy->DeepCopy(inout);
    VTK_CREATE(vtkDelaunay2D, delaunay);
    VTK_CONNECT(delaunay, copy);

    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    vtkPolyData* output=delaunay->GetOutput();

    VTK_CREATE(vtkCellArray, empty_array);
    inout->SetPolys(empty_array);
    inout->GetPolys()->DeepCopy(delaunay->GetOutput()->GetPolys());
#else
    VTK_CREATE(vtkDelaunay2D, delaunay);
    VTK_CONNECT(delaunay, inout);
    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    inout->SetPolys(delaunay->GetOutput()->GetPolys());
#endif

    return inout;
}

inline vtkPolyData* add_mesh2d(vtkPolyData* inout)
{
    VTK_CREATE(vtkPolyData, copy);
    copy->DeepCopy(inout);

    VTK_CREATE(vtkDelaunay2D, delaunay);
    VTK_CONNECT(delaunay, copy);
    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    inout->GetPolys()->DeepCopy(delaunay->GetOutput()->GetPolys());
    return inout;
}

template<typename ForwardContainer>
inline vtkSmartPointer<vtkPolyData>&
add_mesh2d(vtkSmartPointer<vtkPolyData>& inout,
           const ForwardContainer& triangles)
{
    typedef typename ForwardContainer::const_iterator  iterator_type;

    VTK_CREATE(vtkCellArray, cells);
    for (iterator_type it=triangles.begin() ; it!=triangles.end() ;) {
        cells->InsertNextCell(3);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    inout->SetPolys(cells);
    return inout;
}

template<typename ForwardContainer>
inline vtkPolyData* add_mesh2d(vtkPolyData* inout,
                               const ForwardContainer& triangles)
{
    typedef typename ForwardContainer::const_iterator  iterator_type;

    add_vertices(inout);
    VTK_CREATE(vtkCellArray, cells);
    for (iterator_type it=triangles.begin() ; it!=triangles.end() ;) {
        cells->InsertNextCell(3);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    inout->SetPolys(cells);
    return inout;
}

/*
template<typename dataset_type, typename Enable=void>
struct cell_links {};

template<typename dataset_type>
struct cell_links<dataset_type, typename std::enable_if<is_structured<dataset_type>::value>::type> {
    typedef dataset_type dataset_t;

    cell_links(dataset_t* dataset) {
        dataset->GetDimensions(_dims);
    }

    void get_neighbors(std::vector<vtkIdType>& neighbors, unsigned int i) {
        constexpr int all_neighbors[26][3] = {
            {-1, -1, -1}, {-1, 0, -1}, {-1, 1, -1},
            { 0, -1, -1}, { 0, 0, -1}, { 0, 1, -1},
            { 1, -1, -1}, { 1, 0, -1}, { 1, 1, -1},
            {-1, -1,  0}, {-1, 0,  0}, {-1, 1,  0},
            { 0, -1,  0}, { 0, 1,  0},
            { 1, -1,  0}, { 1, 0,  0}, { 1, 1,  0},
            {-1, -1,  1}, {-1, 0,  1}, {-1, 1,  1},
            { 0, -1,  1}, { 0, 0,  1}, { 0, 1,  1},
            { 1, -1,  1}, { 1, 0,  1}, { 1, 1,  1},
        };
        int c[3];
        c[0] = i % _dims[0];
        int j = i / _dims[0];
        c[1] = j % _dims[1];
        c[2] = j / _dims[1];
        for (int n=0; n<26; ++n) {
            const int *ids = all_neighbors[n];
            bool skip = false;
            for (int d=0; d<3 && !skip; ++d) {
                if ( (c[d]==0 && ids[d]==-1) ||
                     (c[d]==_dims[d]-1 && ids[d]==1) )
                    skip = true;
            }
            if (skip) continue;
            int _c[3] = { c[0]+ids[0], c[1]+ids[1], c[2]+ids[2] };
            neighbors.push_back(_c[0] + _dims[0]*(_c[1] + _dims[1]*_c[2]));
        }
    }

    int _dims[3];
};
*/

template<typename DataSet>
struct cell_links {
    typedef DataSet dataset_t;
    typedef cell_links<dataset_t> self_t;

    cell_links(dataset_t* dataset) {
        VTK_INIT(vtkStaticCellLinks, _links);
        _dataset = dataset;
        std::cout << "Building cell links\n";
        _links->Initialize();
        _links->BuildLinks(_dataset);
        std::cout << "Cell links built\n";
        _links->PrintSelf(std::cout, vtkIndent(0));
        // _links = vtkCellLinks::SafeDownCast(dataset->GetCellLinks());
    }

    void get_neighbors(std::vector<vtkIdType>& neighbors, unsigned int i) {
        int nneighs = _links->GetNumberOfCells(i);
        // std::cout << "vertex " << i << " is surrounded by " << nneighs << " cells\n";
        std::set<size_t> point_ids;
        vtkIdType* cell_ids = _links->GetCells(i);
        for (size_t n=0; n<nneighs; ++n) {
            vtkIdType cellid = cell_ids[n];
            VTK_CREATE(vtkIdList, cellpts);
            // std::cout << "from neighbor cell " << n << "(" << cellid << "): ";
            _dataset->GetCellPoints(cellid, cellpts);
            for (size_t k=0; k<cellpts->GetNumberOfIds(); ++k) {
                vtkIdType id = cellpts->GetId(k);
                if (id != i) {
                    point_ids.insert(id);
                    // std::cout << id << ", ";
                }
                else {
                    // std::cout << "(" << i << "), ";
                }
            }
            // std::cout << '\n';
        }
        neighbors = std::vector<vtkIdType>(point_ids.begin(), point_ids.end());
        // std::cout << "there are " << neighbors.size() << " unique vertex neighbors\n";
    }

    VTK_SMART(vtkStaticCellLinks) _links;
    VTK_SMART(dataset_t) _dataset;
};

template<typename dataset_type>
VTK_SMART(dataset_type)
add_derivative_impl(VTK_SMART(dataset_type) inout,
               const std::string& value_name,
               const std::string& deriv_name="") {
    typedef Eigen::MatrixXd matrix_t;
    typedef Eigen::VectorXd vector_t;
    typedef Eigen::Vector<double, 3> point_t;

    /*
        Quadratic fit in 3D

        [1, x, y, z, x^2, xy, xz, y^2, yz, z^2]

        d/dx [0, 1, 0, 0, 2x, y, z, 0, 0, 0]
        d/dy [0, 0, 1, 0, 0, x, 0, 2y, z, 0]
        d/dz [0, 0, 0, 1, 0, 0, x, 0, y, 2z]

        d2/dx2  [0, 0, 0, 0, 2, 0, 0, 0, 0, 0]
        d2/dxdy [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
        d2/dxdz [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
        d2/dy2  [0, 0, 0, 0, 0, 0, 0, 2, 0, 0]
        d2/dydz [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
        d2dz2   [0, 0, 0, 0, 0, 0, 0, 0, 0, 2]

        f(x0+p) - f(x0) = a.x + b.y + c.z + d.x^2 + e.xy + f.xz + g.y^2 + h.yz + i*z^2
        d/dx [] = a + 2dx + ey + fz
        d/dy [] = b + ex + 2gy + hz
        d/dz [] = c + fx + hy + 2iz
        d2/dx2 [] = 2d
        d2/dxdy [] = e
        d2/dxdz [] = f
        d2/dy2 [] = 2g
        d2/dydz [] = h
        d2/dz2 [] = 2i
    */

    cell_links<dataset_type> links(inout);
    size_t npts = inout->GetNumberOfPoints();

    int dummy;
    VTK_SMART(vtkDataArray) array = inout->GetPointData()->GetArray(value_name.c_str(), dummy);
    VTK_CREATE(vtkDoubleArray, derivative);
    if (!deriv_name.empty()) {
        derivative->SetName(deriv_name.c_str());
    }
    else {
        std::string n = "derivative of ";
        n += array->GetName();
        derivative->SetName(n.c_str());
    }
    int ncomp = array->GetNumberOfComponents();
    derivative->SetNumberOfComponents(3*ncomp);
    derivative->SetNumberOfTuples(npts);
    vtkIdType _count = 0;
    vector_t tensor(3*ncomp);

    for (size_t i=0; i<npts; ++i) {
        double tmp[12];
        inout->GetPoint(i, tmp);
        point_t x0(tmp[0], tmp[1], tmp[2]);
        array->GetTuple(i, tmp);
        vector_t f0(ncomp);
        std::copy(tmp, tmp + ncomp, &f0(0));
        std::vector<vtkIdType> point_ids;
        links.get_neighbors(point_ids, i);
        std::vector<point_t> points;
        std::vector<vector_t> values;
        for (auto it=point_ids.begin(); it!=point_ids.end(); ++it) {
            vtkIdType ptid = *it;
            inout->GetPoint(ptid, tmp);
            points.push_back(point_t(tmp[0], tmp[1], tmp[2])-x0);
            array->GetTuple(ptid, tmp);
            vector_t f(ncomp);
            std::copy(tmp, tmp + ncomp, &f(0));
            values.push_back(f-f0);
        }
        // compute weighted linear fit of surrounding values
        matrix_t A(points.size(), 3);
        matrix_t rhs(points.size(), ncomp);
        for (size_t n=0; n<points.size(); ++n) {
            double dsq = points[n].squaredNorm();
            A(n,0) = points[n][0]/(dsq+1.0e-6);
            A(n,1) = points[n][1]/(dsq+1.0e-6);
            A(n,2) = points[n][2]/(dsq+1.0e-6);
            for (int k=0; k<ncomp; ++k) {
                rhs(n,k) = values[n][k]/(dsq+1.0e-6);
            }
        }
        if (points.size() == 0) {
            std::fill(&tensor(0), &tensor(0)+ncomp*3, (double)0);
            // tensor = vector_t::Zero();
            std::cerr << "WARNING: Vertex #" << i << " has no neighbors!";
            double q[3];
            inout->GetPoint(i, q);
            std::cerr << " Position: (" << q[0] << ", " << q[1] << "," << q[2] << ")\n";
        }
        else {
            Eigen::JacobiSVD<matrix_t> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
            matrix_t fit = svd.solve(rhs);
            // the columns of "fit" are the rows of the jacobian
            for (int col=0; col<ncomp; ++col) {
                for (int d=0; d<3; ++d) {
                    tensor[3*col+d] = fit(d, col);
                }
            }
        }

        derivative->SetTypedTuple(_count++, &tensor(0));
    }
    inout->GetPointData()->AddArray(derivative);
    return inout;
}

template<typename dataset_type>
VTK_SMART(dataset_type)
add_derivative(VTK_SMART(dataset_type) inout,
               const std::string& value_name,
               const std::string& deriv_name="") {
    if (vtkUnstructuredGrid::SafeDownCast(inout) != nullptr) {
        VTK_SMART(vtkUnstructuredGrid) ugrid = vtkUnstructuredGrid::SafeDownCast(inout);
        add_derivative_impl(ugrid, value_name, deriv_name);
    }
    else if (vtkRectilinearGrid::SafeDownCast(inout) != nullptr) {
        VTK_SMART(vtkRectilinearGrid) rgrid = vtkRectilinearGrid::SafeDownCast(inout);
        add_derivative_impl(rgrid, value_name, deriv_name);
    }
    else if (vtkStructuredGrid::SafeDownCast(inout) != nullptr) {
        VTK_SMART(vtkStructuredGrid) sgrid = vtkStructuredGrid::SafeDownCast(inout);
        add_derivative_impl(sgrid, value_name, deriv_name);
    }
    else if (vtkImageData::SafeDownCast(inout) != nullptr) {
        VTK_SMART(vtkImageData) igrid = vtkImageData::SafeDownCast(inout);
        add_derivative_impl(igrid, value_name, deriv_name);
    }
    else {
        std::cerr << "Unsupported dataset type!";
    }
    return inout;
}

template<typename dataset_type>
VTK_SMART(dataset_type)
add_jacobian(VTK_SMART(dataset_type) inout, const std::string& name="Jacobian")
{
    if (inout->GetPointData()->GetVectors() == nullptr) {
        std::cerr << "ERROR: No vector field available for Jacobian computation!\n";
        throw std::runtime_error("ERROR: No vector field available for Jacobian computation!\n");
    }
    inout = add_derivative(inout, inout->GetPointData()->GetVectors()->GetName(), name);

    return inout;
}

template<typename dataset_type>
VTK_SMART(dataset_type)
add_gradient(VTK_SMART(dataset_type) inout, const std::string& name="Gradient")
{
    if (inout->GetPointData()->GetScalars() == nullptr) {
        std::cerr << "ERROR: No scalar field available for gradient computation!\n";
        throw std::runtime_error("ERROR: No scalar field available for gradient computation!\n");
    }
    inout = add_derivative(inout, inout->GetPointData()->GetScalars()->GetName(), name);

    return inout;
}

template<typename dataset_type>
VTK_SMART(dataset_type)
add_hessian(VTK_SMART(dataset_type) inout, const std::string& gradient_name="Gradient", const
            std::string& hessian_name="Hessian") {
    if (inout->GetPointData()->GetArray(gradient_name.c_str()) == nullptr) {
        try {
            add_gradient(inout, gradient_name);
        }
        catch (std::exception& e) {
            std::cerr << "Exception caught in add_hessian:\n";
            std::cerr << e.what() << '\n';
            throw;
        }
    }

    inout = add_derivative(inout, "Gradient", "Hessian");
    return inout;
}

inline vtkSmartPointer<vtkUnstructuredGrid>
add_vorticity(vtkSmartPointer<vtkUnstructuredGrid> inout, const std::string& jacobian_name="Jacobian")
{
    typedef typename vtk_array_traits<double>::array_type array_type;
    VTK_SMART(vtkDataArray) jacobian = inout->GetPointData()->GetArray(jacobian_name.c_str());
    if (jacobian.Get() == nullptr) {
        add_jacobian(inout, jacobian_name);
        jacobian = inout->GetPointData()->GetArray(jacobian_name.c_str());
    }

    size_t npts = inout->GetNumberOfPoints();
    VTK_CREATE(array_type, vorticity);
    vorticity->SetNumberOfComponents(3);
    vorticity->SetNumberOfTuples(npts);
    vorticity->SetName("vorticity");
    double vec[3] = {0, 0, 0};
    double tens[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (size_t i=0; i<npts; ++i) {
        jacobian->GetTuple(i, tens);
        vec[0] = tens[7] - tens[5];
        vec[1] = tens[2] - tens[6];
        vec[2] = tens[3] - tens[1];
        vorticity->SetTypedTuple(i, vec);
    }
    inout->GetPointData()->AddArray(vorticity);
    return inout;
}

inline vtkSmartPointer<vtkUnstructuredGrid>
add_lambda2(vtkSmartPointer<vtkUnstructuredGrid> inout, const std::string& jacobian_name="Jacobian")
{
    typedef typename vtk_array_traits<double>::array_type array_type;
    VTK_SMART(vtkDataArray) jacobian = inout->GetPointData()->GetArray(jacobian_name.c_str());
    if (jacobian.Get() == nullptr) {
        add_jacobian(inout, jacobian_name);
        jacobian = inout->GetPointData()->GetArray(jacobian_name.c_str());
    }

    size_t npts = inout->GetNumberOfPoints();
    VTK_CREATE(array_type, lambda2);
    lambda2->SetNumberOfComponents(1);
    lambda2->SetNumberOfTuples(npts);
    lambda2->SetName("lambda2");

    Eigen::Matrix3d J, S, Omega;
    for (size_t i=0; i<npts; ++i) {
        jacobian->GetTuple(i, reinterpret_cast<double*>(&J(0,0)));
        J.transposeInPlace();
        S = 0.5*(J + J.transpose());
        Omega = 0.5*(J - J.transpose());
        S *= S;
        Omega *= Omega;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(S + Omega);
        Eigen::Vector3d evalue = solver.eigenvalues();      // compute eigenvalues

        lambda2->SetTuple(i, &evalue[1]);
    }
    inout->GetPointData()->AddArray(lambda2);
    return inout;
}


// NOTE: Tetrahedras are not supported by vtkPolyData. Hence the output data
//       type has to be vtkUnstructuredGrid, which prevents a direct change
//       of the input parameter unless it is itself ot type vtkUnstructuredGrid
inline vtkSmartPointer<vtkUnstructuredGrid>&
add_mesh3d(vtkSmartPointer<vtkUnstructuredGrid>& inout)
{
    VTK_PTR(vtkDelaunay3D, delaunay);
    VTK_CONNECT(delaunay, inout.GetPointer());
    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    inout->GetCells()->DeepCopy(delaunay->GetOutput()->GetCells());
    return inout;
}

inline vtkUnstructuredGrid* add_mesh3d(vtkUnstructuredGrid* inout)
{
    VTK_PTR(vtkDelaunay3D, delaunay);
    VTK_CONNECT(delaunay, inout);
    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    inout->GetCells()->DeepCopy(delaunay->GetOutput()->GetCells());
    return inout;
}

template<typename ForwardContainer>
inline vtkSmartPointer<vtkUnstructuredGrid>&
add_mesh3d(vtkSmartPointer<vtkUnstructuredGrid>& inout,
           const ForwardContainer& tetrahedra)
{
    typedef typename ForwardContainer::const_iterator  iterator_type;

    VTK_CREATE(vtkCellArray, cells);
    for (iterator_type it=tetrahedra.begin() ; it!=tetrahedra.end() ;) {
        cells->InsertNextCell(4);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    inout->SetCells(VTK_TETRA, cells);
    return inout;
}

template<typename ForwardContainer>
inline vtkUnstructuredGrid* add_mesh3d(vtkUnstructuredGrid* inout,
                                       const ForwardContainer& tetrahedra)
{
    typedef typename ForwardContainer::const_iterator  iterator_type;

    VTK_PTR(vtkCellArray, cells);
    for (iterator_type it=tetrahedra.begin() ; it!=tetrahedra.end() ;) {
        cells->InsertNextCell(4);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    inout->SetCells(VTK_TETRA, cells);
    return inout;
}

inline vtkUnstructuredGrid* create_mesh3d(vtkSmartPointer<vtkPolyData>& _in)
{
    VTK_PTR(vtkDelaunay3D, delaunay);
    VTK_CONNECT(delaunay, _in);
    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    return delaunay->GetOutput();
}

inline vtkUnstructuredGrid* create_mesh3d(vtkPolyData* _in)
{
    add_vertices(_in);
    VTK_PTR(vtkDelaunay3D, delaunay);
    VTK_CONNECT(delaunay, _in);
    delaunay->BoundingTriangulationOff();
    delaunay->Update();
    return delaunay->GetOutput();
}

template<typename ForwardContainer>
inline vtkUnstructuredGrid* create_mesh3d(vtkSmartPointer<vtkPolyData>& _in,
                                          const ForwardContainer& tetrahedra)
{
    typedef typename ForwardContainer::const_iterator  iterator_type;

    VTK_CREATE(vtkCellArray, cells);
    for (iterator_type it=tetrahedra.begin() ; it!=tetrahedra.end() ;) {
        cells->InsertNextCell(4);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    VTK_PTR(vtkUnstructuredGrid, _out);
    _out->GetPoints()->DeepCopy(_in->GetPoints());
    _out->GetPointData()->DeepCopy(_in->GetPointData());
    _out->SetCells(VTK_TETRA, cells);
    return _out;
}

template<typename ForwardContainer>
inline vtkUnstructuredGrid* create_mesh3d(vtkPolyData* _in,
                                          const ForwardContainer& tetrahedra)
{
    typedef typename ForwardContainer::const_iterator  iterator_type;

    VTK_CREATE(vtkCellArray, cells);
    for (iterator_type it=tetrahedra.begin() ; it!=tetrahedra.end() ;) {
        cells->InsertNextCell(4);
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
        cells->InsertCellPoint((*it++));
    }
    VTK_PTR(vtkUnstructuredGrid, _out);
    _out->GetPoints()->DeepCopy(_in->GetPoints());
    _out->GetPointData()->DeepCopy(_in->GetPointData());
    _out->SetCells(VTK_TETRA, cells);
    return _out;
}

template<typename Position_, typename Value_, typename Int_>
inline vtkPolyData* create_mesh(const std::vector<Position_>& pos,
                                const std::string& name,
                                const std::vector<Value_>& scalars,
                                const std::vector<Int_>& triangles)
{
    assert(pos.size() == scalars.size());

    vtkPolyData* _out = make_points(pos);
    add_scalars(_out, scalars);
    _out->GetPointData()->GetScalars()->SetName(name.c_str());
    add_mesh2d(_out);
    return _out;
}

template<typename Position_, typename Vector_, typename Int_,
         typename = typename std::enable_if<spurt::is_array<Vector_>::value>::type >
inline vtkPolyData* create_mesh(const std::vector<Position_>& pos,
                                const std::string& name,
                                const std::vector<Vector_>& vectors,
                                const std::vector<Int_>& triangles)
{
    assert(pos.size() == vectors.size());
    vtkPolyData* _out = make_points(pos);
    add_vectors(_out, vectors);
    _out->GetPointData()->GetScalars()->SetName(name.c_str());
    add_mesh2d(_out, triangles);
    return _out;
}

template<typename T>
inline T* clip_polydata(T* ds, vtkPlane* plane)
{
    VTK_CREATE(vtkClipPolyData, clip);
    VTK_CONNECT(clip, ds);
    clip->SetClipFunction(plane);
    clip->GenerateClipScalarsOff();
    clip->GenerateClippedOutputOn();
    clip->SetValue(0);
    clip->Update();

    T* out = T::New();
    out->DeepCopy(clip->GetOutput());
    return out;
}

inline vtkPlane* make_plane(const nvis::vec3& normal, const nvis::vec3& x)
{
    VTK_PTR(vtkPlane, plane);
    plane->SetOrigin(x[0], x[1], x[2]);
    plane->SetNormal(normal[0], normal[1], normal[2]);
    return plane;
}

inline size_t data_type_size(const vtkDataArray* _in)
{
    int type = const_cast<vtkDataArray*>(_in)->GetDataType();
    switch (type) {
        case VTK_BIT:
        case VTK_CHAR:
        case VTK_SIGNED_CHAR:
        case VTK_UNSIGNED_CHAR:
            return 1;
        case VTK_SHORT:
        case VTK_UNSIGNED_SHORT:
            return sizeof(short);
        case VTK_INT:
        case VTK_UNSIGNED_INT:
            return sizeof(int);
        case VTK_LONG:
        case VTK_UNSIGNED_LONG:
            return sizeof(long);
        case VTK_LONG_LONG:
        case VTK_UNSIGNED_LONG_LONG:
            return sizeof(long long);
        case VTK_FLOAT:
            return sizeof(float);
        case VTK_DOUBLE:
            return sizeof(double);
        default:
            throw std::runtime_error("data type not supported");
    }
}

inline bool is_integral(const vtkDataArray* _in)
{
    int in_type = const_cast<vtkDataArray*>(_in)->GetDataType();
    switch (in_type) {
        case VTK_BIT:
        case VTK_CHAR:
        case VTK_SIGNED_CHAR:
        case VTK_UNSIGNED_CHAR:
        case VTK_SHORT:
        case VTK_UNSIGNED_SHORT:
        case VTK_INT:
        case VTK_UNSIGNED_INT:
        case VTK_LONG:
        case VTK_UNSIGNED_LONG:
        case VTK_ID_TYPE:
        case VTK_LONG_LONG:
        case VTK_UNSIGNED_LONG_LONG:
        // case VTK_INT64:
        // case VTK_UNSIGNED_INT64:
            return true;
        default:
            return false;
    }
}

inline bool is_floating_point(const vtkDataArray* _in)
{
    int in_type = const_cast<vtkDataArray*>(_in)->GetDataType();
    return (in_type == VTK_FLOAT || in_type == VTK_DOUBLE );
}

inline bool is_unsigned(const vtkDataArray* _in)
{
    int in_type = const_cast<vtkDataArray*>(_in)->GetDataType();
    switch (in_type) {
        case VTK_UNSIGNED_CHAR:
        case VTK_CHAR:
        case VTK_UNSIGNED_SHORT:
        case VTK_UNSIGNED_INT:
        case VTK_UNSIGNED_LONG_LONG:
            return true;
        default:
            return false;
    }
}

inline bool is_signed(const vtkDataArray* _in)
{
    return ((is_integral(_in) || is_floating_point(_in)) &&
            !is_unsigned(_in));
}

template<typename TypeOut, typename TypeIn>
inline vtkDataArray* explicit_cast(const vtkDataArray* _in)
{
    typedef TypeIn                                  in_data_type;
    typedef TypeOut                                 out_data_type;
    typedef vtk_array_traits<TypeIn>                in_array_traits;
    typedef vtk_array_traits<TypeOut>               out_array_traits;
    typedef typename in_array_traits::array_type    in_array_type;
    typedef typename out_array_traits::array_type   out_array_type;

    VTK_PTR(out_array_type, _out);
    const in_array_type* cast_in = static_cast<const in_array_type*>(_in);
    size_t tuple_size = const_cast<in_array_type*>(cast_in)->GetNumberOfComponents();
    size_t array_size = const_cast<in_array_type*>(cast_in)->GetNumberOfTuples();
    _out->SetNumberOfComponents(tuple_size);
    _out->SetNumberOfTuples(array_size);
    in_data_type* tuple_in   = new in_data_type[tuple_size];
    out_data_type* tuple_out = new out_data_type[tuple_size];
    for (size_t i=0 ; i<array_size ; ++i) {
        const_cast<in_array_type*>(cast_in)->GetTypedTuple(i, tuple_in);
        for (size_t j=0 ; j<tuple_size ; ++j) {
            tuple_out[j] = static_cast<out_data_type>(tuple_in[j]);
        }
        _out->SetTypedTuple(i, tuple_out);
    }
    delete[] tuple_in;
    delete[] tuple_out;
    return static_cast<vtkDataArray*>(_out);
}

template<typename T>
inline vtkDataArray* implicit_cast(const vtkDataArray* _in)
{
    int input_type = const_cast<vtkDataArray*>(_in)->GetDataType();
    switch (input_type) {
        case VTK_CHAR:
        case VTK_SIGNED_CHAR:
            return explicit_cast<T, char              >(_in);
        case VTK_UNSIGNED_CHAR:
            return explicit_cast<T, unsigned char     >(_in);
        case VTK_SHORT:
            return explicit_cast<T, short             >(_in);
        case VTK_UNSIGNED_SHORT:
            return explicit_cast<T, unsigned short    >(_in);
        case VTK_INT:
            return explicit_cast<T, int               >(_in);
        case VTK_UNSIGNED_INT:
            return explicit_cast<T, unsigned int      >(_in);
        case VTK_LONG:
            return explicit_cast<T, long              >(_in);
        case VTK_UNSIGNED_LONG:
            return explicit_cast<T, unsigned long     >(_in);
        case VTK_LONG_LONG:
            return explicit_cast<T, long long         >(_in);
        case VTK_UNSIGNED_LONG_LONG:
            return explicit_cast<T, unsigned long long>(_in);
        case VTK_FLOAT:
            return explicit_cast<T, float             >(_in);
        case VTK_DOUBLE:
            return explicit_cast<T, double            >(_in);
        default:
            throw std::runtime_error("unsupported input data type");
    }
}

template<typename TypeOut, typename TypeIn>
inline vtkDataArray* explicit_scale(const vtkDataArray* _in,
                                    const std::pair<double, double>& range_in,
                                    const std::pair<double, double>& range_out)
{
    typedef TypeIn                                  in_data_type;
    typedef TypeOut                                 out_data_type;
    typedef vtk_array_traits<TypeIn>                in_array_traits;
    typedef vtk_array_traits<TypeOut>               out_array_traits;
    typedef typename in_array_traits::array_type    in_array_type;
    typedef typename out_array_traits::array_type   out_array_type;

    VTK_PTR(out_array_type, _out);
    const in_array_type* cast_in = static_cast<const in_array_type*>(_in);
    size_t tuple_size = const_cast<in_array_type*>(cast_in)->GetNumberOfComponents();
    size_t array_size = const_cast<in_array_type*>(cast_in)->GetNumberOfTuples();
    _out->SetNumberOfComponents(tuple_size);
    _out->SetNumberOfTuples(array_size);
    in_data_type* tuple_in   = new in_data_type[tuple_size];
    out_data_type* tuple_out = new out_data_type[tuple_size];

    const double& min_in  = range_in.first;
    const double& min_out = range_out.first;
    const double scale = (range_out.second-range_out.first)/(range_in.second-range_in.first);
    for (size_t i=0 ; i<array_size ; ++i) {
        const_cast<in_array_type*>(cast_in)->GetTypedTuple(i, tuple_in);
        for (size_t j=0 ; j<tuple_size ; ++j) {
            double val = tuple_in[j];
            tuple_out[j] = static_cast<out_data_type>(min_out + (val-min_in)*scale);
        }
        _out->SetTypedTuple(i, tuple_out);
    }
    delete[] tuple_in;
    delete[] tuple_out;
    return static_cast<vtkDataArray*>(_out);
}

template<typename T>
inline vtkDataArray* implicit_scale(const vtkDataArray* _in,
                                    const std::pair<double, double>& range_in,
                                    const std::pair<double, double>& range_out)
{
    int input_type = const_cast<vtkDataArray*>(_in)->GetDataType();
    switch (input_type) {
        case VTK_CHAR:
        case VTK_SIGNED_CHAR:
            return explicit_scale<T, char              >(_in, range_in, range_out);
        case VTK_UNSIGNED_CHAR:
            return explicit_scale<T, unsigned char     >(_in, range_in, range_out);
        case VTK_SHORT:
            return explicit_scale<T, short             >(_in, range_in, range_out);
        case VTK_UNSIGNED_SHORT:
            return explicit_scale<T, unsigned short    >(_in, range_in, range_out);
        case VTK_INT:
            return explicit_scale<T, int               >(_in, range_in, range_out);
        case VTK_UNSIGNED_INT:
            return explicit_scale<T, unsigned int      >(_in, range_in, range_out);
        case VTK_LONG:
            return explicit_scale<T, long              >(_in, range_in, range_out);
        case VTK_UNSIGNED_LONG:
            return explicit_scale<T, unsigned long     >(_in, range_in, range_out);
        case VTK_LONG_LONG:
            return explicit_scale<T, long long         >(_in, range_in, range_out);
        case VTK_UNSIGNED_LONG_LONG:
            return explicit_scale<T, unsigned long long>(_in, range_in, range_out);
        case VTK_FLOAT:
            return explicit_scale<T, float             >(_in, range_in, range_out);
        case VTK_DOUBLE:
            return explicit_scale<T, double            >(_in, range_in, range_out);
        default:
            throw std::runtime_error("unsupported input data type");
    }
}

template<typename T>
inline vtkDataArray* convert(const vtkDataArray* _in)
{
    // Here is how we handle the various cases that may present themselves
    // Note: cases are evaluated in this order
    // 1) integral to floating point         -> cast
    // 2) floating point to floating point   -> cast
    // 3) floating point to integral         -> quantize
    // 4) increase precision (both integral) -> shift
    // 5) decrease precision (both integral) -> quantize

    bool input_is_integral        = is_integral(_in);
    bool output_is_integral       = std::is_integral<T>::value;
    bool input_is_floating_point  = is_floating_point(_in);
	bool output_is_floating_point = std::is_floating_point<T>::value;

    if (!input_is_integral && !input_is_floating_point) {
        throw std::runtime_error("data type of input array is not arithmetic");
    }
    if (!output_is_integral && !output_is_floating_point) {
        throw std::runtime_error("output data type is not arithmetic");
    }

    // 1) integral to floating point         -> cast
    if (input_is_integral && output_is_floating_point) {
        return implicit_cast<T>(_in);
    }
    // 2) floating point to floating point   -> cast
    else if (input_is_floating_point && output_is_floating_point) {
        return implicit_cast<T>(_in);
    }

    double* range = const_cast<vtkDataArray*>(_in)->GetRange();
    std::pair<double, double> range_in, range_out;
    range_in.first = range[0];
    range_in.second = range[1];
    range_out.first = std::numeric_limits<T>::min();
    range_out.second = std::numeric_limits<T>::max();

    // 3) floating point to integral         -> quantize
    if (input_is_floating_point && output_is_integral) {
        return implicit_scale<T>(_in, range_in, range_out);
    }

    // 4a) Input range fits in output range  -> cast
    if (range_in.first >= range_out.first &&
            range_in.second <= range_out.second) {
        implicit_cast<T>(_in);
    }
    // 4b) Input fits after translation      -> shift
    else if (range_in.second-range_in.first <= range_out.second-range_out.first) {
        std::pair<double, double> partial_range_out=std::make_pair(range_out.first, range_out.first+(range_in.second-range_in.first));
        return implicit_scale<T>(_in, range_in, partial_range_out);
    }
    // 5) Output range does not fit          -> quantize
    else {
        return implicit_scale<T>(_in, range_in, range_out);
    }

    return 0; // superfluous fallback options for picky compilers
}

template<typename T>
inline vtkDataArray* merge(const vtkDataArray* _in1,
                           const vtkDataArray* _in2,
                           const vtkDataArray* _in3 = 0,
                           const vtkDataArray* _in4 = 0,
                           const vtkDataArray* _in5 = 0,
                           const vtkDataArray* _in6 = 0,
                           const vtkDataArray* _in7 = 0,
                           const vtkDataArray* _in8 = 0,
                           const vtkDataArray* _in9 = 0)
{
    typedef T                                   data_type;
    typedef vtk_array_traits<T>                 array_traits;
    typedef typename array_traits::array_type   array_type;

    std::vector<const vtkDataArray* > to_merge;
    to_merge.push_back(_in1);
    to_merge.push_back(_in2);
    if (_in3) {
        to_merge.push_back(_in3);
    }
    if (_in4) {
        to_merge.push_back(_in4);
    }
    if (_in5) {
        to_merge.push_back(_in5);
    }
    if (_in6) {
        to_merge.push_back(_in6);
    }
    if (_in7) {
        to_merge.push_back(_in7);
    }
    if (_in8) {
        to_merge.push_back(_in8);
    }
    if (_in9) {
        to_merge.push_back(_in9);
    }

    // ensure consistency across arrays
    size_t size          = const_cast<vtkDataArray*>(_in1)->GetNumberOfTuples();
    size_t nb_components = 0;
    size_t nb_arrays     = to_merge.size();

    std::cout << "there are " << nb_arrays << " arrays in input\n";

    std::vector<array_type*> converted(nb_arrays);
    for (size_t i=0 ; i<nb_arrays ; ++i) {
        if (const_cast<vtkDataArray*>(to_merge[i])->GetNumberOfTuples() !=
            size) {
            throw std::runtime_error("array sizes do not match");
        }
        converted[i] = array_type::SafeDownCast(convert<T>(to_merge[i]));
        nb_components += const_cast<vtkDataArray*>(to_merge[i])->GetNumberOfComponents();
    }

    std::cout << "total number of components is " << nb_components << '\n';

    VTK_CREATE(array_type, array);
    array->SetNumberOfTuples(size);
    array->SetNumberOfComponents(nb_components);
    data_type* tuple_out = new data_type[nb_components];
    data_type* tuple_in  = new data_type[nb_components];
    for (size_t i=0 ; i<size ; ++i) {
        size_t j=0;
        for (size_t n=0 ; n<nb_arrays ; ++n) {
            converted[n]->GetTypedTuple(i, tuple_in);
            for (size_t k=0 ; k<converted[n]->GetNumberOfComponents() ; ++k) {
                tuple_out[j++] = tuple_in[k];
            }
        }
        array->SetTypedTuple(i, tuple_out);
    }
    return static_cast<vtkDataArray*>(array);
}




} // namespace vtk_utils


#endif
