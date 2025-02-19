#ifndef __GARCIA_VIS_HELPER_HPP__
#define __GARCIA_VIS_HELPER_HPP__

#include <vector>
#include <map>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <ostream>
#include <list>
#include "vtk_utils.hpp"

#include "vtkPolyDataReader.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkPolyDataNormals.h"
#include "vtkCommand.h"
#include "vtkInteractorObserver.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkCylinderSource.h"
#include "vtkWindowToImageFilter.h"
#include "vtkTIFFWriter.h"
#include "vtkDataSetReader.h"
#include "vtkTubeFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkCubeSource.h"
#include "vtkContourFilter.h"
#include "vtkStructuredPointsReader.h"
#include "vtkScalarBarActor.h"
#include "vtkColorTransferFunction.h"
#include "vtkTextProperty.h"
#include "vtkStructuredPoints.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"

namespace Garcia_vis_helper {

struct dataset_info {
    std::string base_dir;
    std::string microstruct;
    std::string mesh_base;
    std::string movie_path;
    std::string stress;
    std::string stress_nrrd;
    std::string id_to_tags;
    std::string stress_span;
    std::string stat_base;
    std::string orientation;
    std::string dfield;
    std::string efield;
    std::string dfield_norm;
    std::string efield_norm;
    std::string dfield_span;
    std::string efield_span;
    nvis::vec3 up, position, focal_point;
    double near, far;
    nvis::vec3 step;
    nvis::bbox3 valid_bounds;
};

// point to grain relationship

typedef std::set<int>   tag_type;

template<typename T>
inline bool included(const tag_type& tags, const T values, int size)
{
    for (int i = 0 ; i < size ; ++i) {
        if (tags.find(values[i]) == tags.end()) {
            return false;
        }
    }
    return true;
}

// display helper

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
    os << "[ ";
    for (typename std::set<T>::const_iterator i = s.begin() ; i != s.end() ; ++i) {
        os << *i << " ";
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& s)
{
    os << "[ ";
    for (typename std::vector<T>::const_iterator i = s.begin() ; i != s.end() ; ++i) {
        os << *i << " ";
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::list<T>& s)
{
    os << "[ ";
    for (typename std::list<T>::const_iterator i = s.begin() ; i != s.end() ; ++i) {
        os << *i << " ";
    }
    os << "]";
    return os;
}

// color related stuff
typedef unsigned char                   uchar;
typedef nvis::fixed_vector<uchar, 3>    uchar_color_type;

template<typename T>
inline T cubic_spline(const T& a, const T& b, double t)
{
    return a + (b - a)*t*t*(3 - 2*t);
}

inline uchar_color_type to_color(const nvis::vec3& v)
{
    return uchar_color_type(255.*v);
}

inline nvis::vec3 updown(double z, double gamma)
{
    nvis::vec3 pole = (z < 0 ? nvis::vec3(0, 0, 1) : nvis::vec3(1, 1, 0));
    double t = pow(fabs(z), gamma);
    return cubic_spline(nvis::vec3(1, 1, 1), pole, t);
}

inline nvis::vec3 vec2col(const nvis::vec3& dir, bool up, double gamma)
{
    nvis::vec3 e(sin(dir[2])*sin(dir[1]),
                 cos(dir[2])*sin(dir[1]),
                 cos(dir[1]));
    double z = e[2];
    
    if (up) {
        return updown(z, gamma);
    } else {
        const nvis::vec3 hue[] = {
            nvis::vec3(0, 0, 1), // 0: blue
            nvis::vec3(0, 1, 1), // 1: cyan
            nvis::vec3(0, 1, 0), // 2: green
            nvis::vec3(0.5, 1, 0),
            nvis::vec3(1, 1, 0), // 2: yellow
            nvis::vec3(1, 0.5, 0), // 3: red
            nvis::vec3(1, 0, 0), // 3: red
            nvis::vec3(1, 0, 1), // 3: red
            nvis::vec3(0, 0, 1)  // 4: blue to close the loop
        };
        
        const double delta = M_PI / 4.;
        const nvis::vec3 top = nvis::vec3(1, 1, 1);
        const nvis::vec3 bottom = nvis::vec3(0.3, 0.3, 0.3);
        
        nvis::vec2 xy(e[0], e[1]);
        double norm = nvis::norm(xy);
        if (norm == 0) {
            return (z < 0 ? bottom : top);
        }
        xy /= norm;
        double alpha = acos(xy[0]);
        if (xy[1] < 0) {
            alpha = 2 * M_PI - alpha;
        }
        int i = floor(alpha / delta);
        double u = alpha / delta - (double)i;
        // nvis::vec3 h = cubic_spline(hue[i], hue[i+1], u);
        nvis::vec3 h = (1. - u) * hue[i] + u * hue[i+1];
        
        double t = pow(fabs(z), gamma);
        nvis::vec3 pole = (z < 0 ? bottom : top);
        
        return cubic_spline(h, pole, t);
    }
}

inline nvis::vec3 brighten(const nvis::vec3& c)
{
    double max = *std::max_element(c.begin(), c.end());
    return c / max;
}

inline nvis::vec3 upcol(const nvis::vec3& dir, double brightness, double gamma, bool transform = false)
{
    nvis::vec3 d = dir / nvis::norm(dir);
    nvis::vec3 e;
    if (transform) e = nvis::vec3(sin(d[2]) * sin(d[1]),
                                      cos(d[2]) * sin(d[1]),
                                      cos(d[1]));
    else {
        e = d;
    }
    double z = e[2];
    const nvis::vec3 colors[] = { nvis::vec3(1, 0, 0),
                                  nvis::vec3(1, 0, 1),
                                  nvis::vec3(0, 0, 1),
                                  nvis::vec3(0, 1, 1),
                                  nvis::vec3(0, 1, 0),
                                  nvis::vec3(0.5, 1, 0),
                                  nvis::vec3(1, 1, 0),
                                  nvis::vec3(1, 0.5, 0),
                                  nvis::vec3(1, 0, 0)
                                };
    const double delta = M_PI / 4.;
    const nvis::vec3 top = nvis::vec3(1, 1, 1);
    
    nvis::vec2 xy(e[0], e[1]);
    double norm = nvis::norm(xy);
    if (norm == 0) {
        return top;
    }
    xy /= norm;
    double alpha = acos(xy[0]);
    if (xy[1] < 0) {
        alpha = 2 * M_PI - alpha;
    }
    
    int sector = floor(alpha / delta);
    if (sector > 7) {
        sector = 0;
    }
    
    double u = alpha / delta - sector;
    nvis::vec3 base_color = ((1. - u) * colors[sector] + u * colors[sector+1]);
    
    double t = pow(fabs(z), gamma);
    
    return ((1 - fabs(z))*brightness + fabs(z)) * brighten(cubic_spline(base_color, top, t));
}

inline nvis::vec3 xyz2rgb(const nvis::vec3& dir, double brightness, double gamma, bool transform = false)
{
    nvis::vec3 d = dir / nvis::norm(dir);
    nvis::vec3 e;
    if (transform)
        e = nvis::vec3(sin(d[2]) * sin(d[1]),
                       cos(d[2]) * sin(d[1]),
                       cos(d[1]));
    else {
        e = d;
    }
    
    nvis::vec3 rgb(fabs(e[2]), fabs(e[0]), fabs(e[1]));
    
    double w = 0.6 + 0.4*pow(fabs(e[2]), gamma);
    nvis::vec3 col = (1-w)*nvis::vec3(1,1,1) + w*rgb;
    col /= nvis::norm(col);
    
    return col;
}

template<typename T>
struct spherical_color_map {
    typedef nvis::fvec3                     color_type;
    typedef T                               data_type;
    typedef nvis::fixed_vector<T, 3>        vec_type;
    
    enum map_kind {
        RED_GREEN_SATURATION,
        RED_GREEN_BLUE,
        GREEN_RED_SATURATION,
        GREEN_RED_BLUE,
        BLUE_YELLOW_SATURATION,
        YELLOW_BLUE_SATURATION,
    };
    
    spherical_color_map(int kind = 0, T gamma = 1) : _kind(kind), _gamma(gamma) {}
    
    color_type operator()(const vec_type& v) const {
        vec_type w = v / nvis::norm(v);
        data_type z = w[2];
        // std::cerr << "w = " << w << ", z = " << z;
        
        vec_type pole, equator;
        switch (_kind) {
            case RED_GREEN_SATURATION:
            case RED_GREEN_BLUE:
                pole = z < 0 ? color_type(0, 1, 0) : color_type(1, 0, 0);
                break;
            case GREEN_RED_SATURATION:
            case GREEN_RED_BLUE:
                pole = z < 0 ? color_type(1, 0, 0) : color_type(0, 1, 0);
                break;
            case BLUE_YELLOW_SATURATION:
                pole = z < 0 ? color_type(1, 1, 0) : color_type(0, 0, 1);
                break;
            case YELLOW_BLUE_SATURATION:
                pole = z < 0 ? color_type(0, 0, 1) : color_type(1, 1, 0);
                break;
            default:
                std::cerr << "invalid color map kind\n";
                exit(-1);
        }
        
        switch (_kind) {
            case RED_GREEN_SATURATION:
            case GREEN_RED_SATURATION:
            case BLUE_YELLOW_SATURATION:
            case YELLOW_BLUE_SATURATION:
                equator = color_type(1, 1, 1);
                break;
            case GREEN_RED_BLUE:
            case RED_GREEN_BLUE:
                equator = color_type(0, 0, 1);
                break;
            default:
                std::cerr << "invalid color map kind\n";
                exit(-1);
        }
        
        double t = pow(fabs(z), _gamma);
        color_type c = cubic_spline(equator, pole, t);
        
        // std::cerr << ", t = " << t << ", color = " << c << '\n';
        return c;
    }
    
    int         _kind;
    data_type   _gamma;
};

// blue - white - red color map for stress field
template<typename T>
struct color_map {
    typedef nvis::fvec3     color_type;
    typedef T               data_type;
    
    enum map_kind {
        DOUBLE_ENDED,
        REDUNDANT_RAINBOW,
        BLUE_TO_YELLOW,
        GREEN_TO_RED,
        BLACK_TO_WHITE,
    };
    
    color_map(const std::vector<T>& values, data_type gamma = 1, bool adaptive = false)
        : __vals(values), __gamma(gamma), __adaptive(adaptive) {
        std::sort(__vals.begin(), __vals.end());
        __min = __vals.front();
        __max = __vals.back();
        
        __width = std::max(fabs(__min), fabs(__max));
    }
    
    color_map(data_type min, data_type max, data_type gamma = 1) : __gamma(gamma) {
        __width = std::max(fabs(min), fabs(max));
    }
    
    color_type double_ended(data_type v) const {
        color_type hue;
        if (v < 0) {
            hue = color_type(0, 0, 1);
        } else {
            hue = color_type(1, 0, 0);
        }
        
        data_type u = fabs(v) / __width;
        data_type x = pow(u, __gamma);
        return (1. - x)*color_type(1, 1, 1) + x*hue;
    }
    
    color_type rainbow(data_type u) const {
        static nvis::vec3 __colors[] = {
            nvis::vec3(0, 0, 0.5),  // 0:  dark blue
            nvis::vec3(0, 0, 1),    // 1:  blue
            nvis::vec3(0, 0.5, 1),  // 2:  sky blue
            nvis::vec3(0, 1, 1),    // 3:  cyan
            nvis::vec3(0, 1, 0.5),  // 4:
            nvis::vec3(0, 1, 0),    // 5:  green
            nvis::vec3(0.5, 1, 0),  // 6:
            nvis::vec3(1, 1, 0),    // 7:  yellow
            nvis::vec3(1, 0.5, 0),  // 8:  orange
            nvis::vec3(1, 0, 0),    // 9:  red
            nvis::vec3(1, 0, 0.5),  // 10:
            nvis::vec3(1, 0, 1),    // 11: magenta
            nvis::vec3(1, 0.5, 1),  // 12:
            nvis::vec3(1, 1, 1)     // 13: white
        };
        
        data_type x = u * (data_type)13;
        int i = floor(x);
        if (i == 13) {
            return __colors[13];
        }
        data_type t = x - (data_type)i;
        return (1. - t)*__colors[i] + t*__colors[i+1];
    }
    
    color_type redundant(data_type v) const {
        data_type u;
        if (__adaptive) {
            int rank = std::distance(__vals.begin(),
                                     std::lower_bound(__vals.begin(), __vals.end(), v));
            u = (data_type)rank / (data_type)__vals.size();
        } else {
            u = (v - __min) / (__max - __min);
        }
        
        data_type value = 0.5 * (1 + pow(u, __gamma));
        color_type cl = rainbow(u);
        return value * cl;
    }
    
    color_type blue_yellow(data_type v) const {
        data_type u;
        if (__adaptive) {
            int rank = std::distance(__vals.begin(),
                                     std::lower_bound(__vals.begin(), __vals.end(), v));
            u = (data_type)rank / (data_type)__vals.size();
        } else {
            u = (v - __min) / (__max - __min);
        }
        
        return (1. - u)*color_type(0, 0, 1) + u*color_type(1, 1, 0);
    }
    
    color_type green_red(data_type v) const {
        data_type u;
        if (__adaptive) {
            int rank = std::distance(__vals.begin(),
                                     std::lower_bound(__vals.begin(), __vals.end(), v));
            u = (data_type)rank / (data_type)__vals.size();
        } else {
            u = (v - __min) / (__max - __min);
        }
        
        return (1. - u)*color_type(0, 1, 0) + u*color_type(1, 0, 0);
    }
    
    color_type black_white(data_type v) const {
        data_type u;
        if (__adaptive) {
            int rank = std::distance(__vals.begin(),
                                     std::lower_bound(__vals.begin(), __vals.end(), v));
            u = (data_type)rank / (data_type)__vals.size();
        } else {
            u = (v - __min) / (__max - __min);
        }
        
        return u*color_type(1, 1, 1);
    }
    
    color_type operator()(data_type v, map_kind k = DOUBLE_ENDED) const {
        switch (k) {
            case DOUBLE_ENDED:
                return double_ended(v);
            case REDUNDANT_RAINBOW:
                return redundant(v);
            case BLUE_TO_YELLOW:
                return blue_yellow(v);
            case GREEN_TO_RED:
                return green_red(v);
            case BLACK_TO_WHITE:
                return black_white(v);
            default:
                throw;
        }
    }
    
    const data_type& width() const {
        return __width;
    }
    
    std::map<data_type, color_type> __map;
    std::vector<data_type>          __vals;
    data_type                       __width, __gamma, __min, __max;
    bool                            __adaptive;
};

vtkSmartPointer<vtkActor> draw_frame(const nvis::bbox3& frame, double radius = 0)
{
    nvis::vec3 diag = frame.size();
    std::cerr << "frame = " << frame << '\n';
    nvis::vec3 logo[] = {
        nvis::vec3(0, 0, 0), nvis::vec3(diag[0], 0, 0),
        nvis::vec3(0, diag[1], 0), nvis::vec3(-diag[0], 0, 0),
        nvis::vec3(0, -diag[1], diag[2]), nvis::vec3(diag[0], 0, 0),
        nvis::vec3(0, diag[1], 0), nvis::vec3(-diag[0], 0, 0)
    };
    vtkIdType __edges[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    
    // create cube outline
    vtkSmartPointer<vtkDoubleArray> __coords = vtkSmartPointer<vtkDoubleArray>::New();
    __coords->SetNumberOfComponents(3);
    nvis::vec3 p = frame.min();
    __coords->InsertNextTuple(p.begin());
    for (int i = 1 ; i < 8 ; ++i) {
        p += logo[i];
        __coords->InsertNextTuple(p.begin());
    }
    vtkSmartPointer<vtkPoints> cube_pts = vtkSmartPointer<vtkPoints>::New();
    cube_pts->SetData(__coords);
    
    vtkSmartPointer<vtkCellArray> cube_lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType n = 2;
    for (int i = 0 ; i < 12 ; ++i) {
        cube_lines->InsertNextCell(n, __edges[i]);
    }
    
    vtkSmartPointer<vtkPolyData> cube_edges = vtkSmartPointer<vtkPolyData>::New();
    cube_edges->SetPoints(cube_pts);
    cube_edges->SetLines(cube_lines);
    
    vtkSmartPointer<vtkTubeFilter> edge_tubes = vtkSmartPointer<vtkTubeFilter>::New();
    
    double r = (radius == 0 ? 0.005* *std::max_element(diag.begin(), diag.end()) : radius);
    edge_tubes->SetInputData(cube_edges);
    edge_tubes->SetRadius(r);
    edge_tubes->SetNumberOfSides(6);
    
    std::cerr << "cube created\n";
    
    vtkSmartPointer<vtkPolyDataMapper> cube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cube_mapper->SetInputConnection(edge_tubes->GetOutputPort());
    vtkSmartPointer<vtkActor> cube_actor = vtkSmartPointer<vtkActor>::New();
    cube_actor->SetMapper(cube_mapper);
    cube_actor->GetProperty()->SetColor(1, 1, 1);
    cube_actor->GetProperty()->SetOpacity(1);
    
    return cube_actor;
}

// hard-coded settings for available datasets
static std::string data_root = "/scratch4/data/Garcia/";

static dataset_info textured_info, untextured_info, mc_info_00, mc_info_06, mc_info_10_09, mc_info_10_05;

static void set_paths()
{
    textured_info.base_dir = data_root + "microstructure/";
    textured_info.mesh_base = textured_info.base_dir + "TexBNKT_190_190_37";
    textured_info.microstruct = textured_info.base_dir + "TexBNKT_190_190_37.nrrd";
    textured_info.stress = textured_info.base_dir + "trace-Estress_TexturedBNKT_b09.vtk";
    textured_info.stress_nrrd = textured_info.base_dir + "trace-Estress_TexturedBNKT_b09.nrrd";
    textured_info.id_to_tags = textured_info.base_dir + "Textured_b09-ids_to_tags.nrrd";
    textured_info.stress_span = textured_info.base_dir + "Textured_b09-grain-span.nrrd";
    textured_info.stat_base = textured_info.base_dir + "Textured_b09";
    textured_info.orientation = textured_info.mesh_base + "-euler.nrrd";
    textured_info.movie_path = textured_info.base_dir + "tex/";
    textured_info.dfield = textured_info.base_dir + "xmt_dfield.vtk";
    textured_info.efield = textured_info.base_dir + "xmt_efield.vtk";
    textured_info.dfield_norm = textured_info.base_dir + "xmt_norm_dfield.vtk";
    textured_info.efield_norm = textured_info.base_dir + "xmt_norm_efield.vtk";
    textured_info.dfield_span = textured_info.base_dir + "xmt_span_dfield.nrrd";
    textured_info.efield_span = textured_info.base_dir + "xmt_span_efield.nrrd";
    textured_info.position = nvis::vec3(-38.2382, -20.2557, 27.1584);
    textured_info.focal_point = nvis::vec3(17.1675, 19.0469, 3.32833);
    textured_info.up = nvis::vec3(0.254529, 0.213122, 0.943289);
    textured_info.near = 17.6447;
    textured_info.far = 140.982;
    textured_info.step = nvis::vec3(0.2, 0.2, 0.4);
    textured_info.valid_bounds.min() = nvis::vec3(10., 10., 5.);
    textured_info.valid_bounds.max() = nvis::vec3(179., 179., 31.);
    
    untextured_info.base_dir = data_root + "microstructure/";
    untextured_info.mesh_base = untextured_info.base_dir + "UnTexBNKT_01_140_140_70";
    untextured_info.microstruct = untextured_info.base_dir + "UnTexBNKT_01_140_140_70.nrrd";
    untextured_info.stress = textured_info.base_dir + "trace-Estress_TexturedBNKT_b09.vtk";
    untextured_info.stress_nrrd = untextured_info.base_dir + "trace-Estress_UntexturedBNKT_b09.nrrd";
    untextured_info.id_to_tags = untextured_info.base_dir + "Untextured_b09-ids_to_tags.nrrd";
    untextured_info.stress_span = untextured_info.base_dir + "Untextured_b09-grain-span.nrrd";
    untextured_info.stat_base = untextured_info.base_dir + "Untextured_b09";
    untextured_info.orientation = untextured_info.mesh_base + "-euler.nrrd";
    untextured_info.movie_path = untextured_info.base_dir + "untex/";
    untextured_info.dfield = untextured_info.base_dir + "xmt_dfield.nrrd";
    untextured_info.efield = untextured_info.base_dir + "xmt_efield.nrrd";
    untextured_info.dfield_norm = untextured_info.base_dir + "xmt_norm_dfield.nrrd";
    untextured_info.efield_norm = untextured_info.base_dir + "xmt_norm_efield.nrrd";
    untextured_info.dfield_span = untextured_info.base_dir + "xmt_span_dfield.nrrd";
    untextured_info.efield_span = untextured_info.base_dir + "xmt_span_efield.nrrd";
    untextured_info.position = nvis::vec3(-11.7274, -4.0038, 11.1954);
    untextured_info.focal_point = nvis::vec3(15.4841, 11.0543, -3.40036);
    untextured_info.up = nvis::vec3(0.352451, 0.239877, 0.904565);
    untextured_info.near = 5.28852;
    untextured_info.far = 39.2914;
    untextured_info.step = nvis::vec3(0.07, 0.07, 0.1);
    untextured_info.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    untextured_info.valid_bounds.max() = nvis::vec3(129., 129., 59.);
    
    mc_info_00.base_dir = data_root + "r00b09/";
    mc_info_00.mesh_base = data_root + "microstructure/MC_588grains_110cubed";
    mc_info_00.microstruct = mc_info_00.mesh_base + ".nrrd";
    mc_info_00.stress = mc_info_00.base_dir + "trace-Estress_CG_588Grains_r00b09.vtk";
    mc_info_00.stress_nrrd = mc_info_00.base_dir + "trace-Estress_CG_588Grains_r00b09.nrrd";
    mc_info_00.id_to_tags = mc_info_00.base_dir + "CG_588Grains_r00b09-ids_to_tags.nrrd";
    mc_info_00.stress_span = mc_info_00.base_dir + "CG_588Grains_r00b09-grain-span.nrrd";
    mc_info_00.stat_base = mc_info_00.base_dir + "CG_588Grains_r00b09";
    mc_info_00.orientation = mc_info_00.base_dir + "orientation.nrrd";
    mc_info_00.movie_path = mc_info_00.base_dir + "movies/";
    mc_info_00.dfield = mc_info_00.base_dir + "xmt_dfield.nrrd";
    mc_info_00.efield = mc_info_00.base_dir + "xmt_efield.nrrd";
    mc_info_00.dfield_norm = mc_info_00.base_dir + "xmt_norm_dfield.nrrd";
    mc_info_00.efield_norm = mc_info_00.base_dir + "xmt_norm_efield.nrrd";
    mc_info_00.dfield_span = mc_info_00.base_dir + "xmt_span_dfield.nrrd";
    mc_info_00.efield_span = mc_info_00.base_dir + "xmt_span_efield.nrrd";
    mc_info_00.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_00.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_00.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_00.near = 157.65;
    mc_info_00.far = 961.655;
    mc_info_00.step = nvis::vec3(2, 2, 2);
    mc_info_00.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_00.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    mc_info_06.base_dir = data_root + "r06b09/";
    mc_info_06.mesh_base = data_root + "microstructure/MC_588grains_110cubed";
    mc_info_06.microstruct = mc_info_06.mesh_base + ".nrrd";
    mc_info_06.stress = mc_info_06.base_dir + "trace-Estress_CG_588Grains_r06b09.vtk";
    mc_info_06.stress_nrrd = mc_info_06.base_dir + "trace-Estress_CG_588Grains_r06b09.nrrd";
    mc_info_06.id_to_tags = mc_info_06.base_dir + "CG_588Grains_r06b09-ids_to_tags.nrrd";
    mc_info_06.stress_span = mc_info_06.base_dir + "CG_588Grains_r06b09-grain-span.nrrd";
    mc_info_06.stat_base = mc_info_06.base_dir + "CG_588Grains_r06b09";
    mc_info_06.orientation = mc_info_06.base_dir + "orientation.nrrd";
    mc_info_06.movie_path = mc_info_06.base_dir + "movies/";
    mc_info_06.dfield = mc_info_06.base_dir + "xmt_dfield.nrrd";
    mc_info_06.efield = mc_info_06.base_dir + "xmt_efield.nrrd";
    mc_info_06.dfield_norm = mc_info_06.base_dir + "xmt_norm_dfield.nrrd";
    mc_info_06.efield_norm = mc_info_06.base_dir + "xmt_norm_efield.nrrd";
    mc_info_06.dfield_span = mc_info_06.base_dir + "xmt_span_dfield.nrrd";
    mc_info_06.efield_span = mc_info_06.base_dir + "xmt_span_efield.nrrd";
    mc_info_06.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_06.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_06.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_06.near = 157.65;
    mc_info_06.far = 961.655;
    mc_info_06.step = nvis::vec3(2, 2, 2);
    mc_info_06.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_06.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    mc_info_10_09.base_dir = data_root + "r10b09/";
    mc_info_10_09.mesh_base = data_root + "microstructure/MC_588grains_110cubed";
    mc_info_10_09.microstruct = mc_info_10_09.mesh_base + ".nrrd";
    mc_info_10_09.stress = mc_info_10_09.base_dir + "trace-Estress_CG_588Grains_r10b09.vtk";
    mc_info_10_09.stress_nrrd = mc_info_10_09.base_dir + "trace-Estress_CG_588Grains_r10b09.nrrd";
    mc_info_10_09.id_to_tags = mc_info_10_09.base_dir + "CG_588Grains_r10b09-ids_to_tags.nrrd";
    mc_info_10_09.stress_span = mc_info_10_09.base_dir + "CG_588Grains_r10b09-grain-span.nrrd";
    mc_info_10_09.stat_base = mc_info_10_09.base_dir + "CG_588Grains_r10b09";
    mc_info_10_09.orientation = mc_info_10_09.base_dir + "orientation.nrrd";
    mc_info_10_09.movie_path = mc_info_10_09.base_dir + "movies/";
    mc_info_10_09.dfield = mc_info_10_09.base_dir + "xmt_dfield.nrrd";
    mc_info_10_09.efield = mc_info_10_09.base_dir + "xmt_efield.nrrd";
    mc_info_10_09.dfield_norm = mc_info_10_09.base_dir + "xmt_norm_dfield.nrrd";
    mc_info_10_09.efield_norm = mc_info_10_09.base_dir + "xmt_norm_efield.nrrd";
    mc_info_10_09.dfield_span = mc_info_10_09.base_dir + "xmt_span_dfield.nrrd";
    mc_info_10_09.efield_span = mc_info_10_09.base_dir + "xmt_span_efield.nrrd";
    mc_info_10_09.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_10_09.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_10_09.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_10_09.near = 157.65;
    mc_info_10_09.far = 961.655;
    mc_info_10_09.step = nvis::vec3(2, 2, 2);
    mc_info_10_09.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_10_09.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    mc_info_10_05.base_dir = data_root + "r10b05/";
    mc_info_10_05.mesh_base = data_root + "microstructure/MC_588grains_110cubed";
    mc_info_10_05.microstruct = mc_info_10_05.mesh_base + ".nrrd";
    mc_info_10_05.stress = mc_info_10_05.base_dir + "trace-Estress_CG_588Grains_r10b05.vtk";
    mc_info_10_05.stress_nrrd = mc_info_10_05.base_dir + "trace-Estress_CG_588Grains_r10b05.nrrd";
    mc_info_10_05.id_to_tags = mc_info_00.base_dir + "CG_588Grains_r00b09-ids_to_tags.nrrd";
    mc_info_10_05.stress_span = mc_info_10_05.base_dir + "CG_588Grains_r10b05-grain-span.nrrd";
    // HACK!!
    mc_info_10_05.stat_base = mc_info_00.base_dir + "CG_588Grains_r00b09";
    // HACK
    mc_info_10_05.orientation = mc_info_10_05.base_dir + "orientation.nrrd";
    mc_info_10_05.movie_path = mc_info_10_05.base_dir + "movies/";
    mc_info_10_05.dfield = mc_info_10_05.base_dir + "xmt_dfield.nrrd";
    mc_info_10_05.efield = mc_info_10_05.base_dir + "xmt_efield.nrrd";
    mc_info_10_05.dfield_norm = mc_info_10_05.base_dir + "xmt_norm_dfield.nrrd";
    mc_info_10_05.efield_norm = mc_info_10_05.base_dir + "xmt_norm_efield.nrrd";
    mc_info_10_05.dfield_span = mc_info_10_05.base_dir + "xmt_span_dfield.nrrd";
    mc_info_10_05.efield_span = mc_info_10_05.base_dir + "xmt_span_efield.nrrd";
    mc_info_10_05.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_10_05.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_10_05.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_10_05.near = 157.65;
    mc_info_10_05.far = 961.655;
    mc_info_10_05.step = nvis::vec3(2, 2, 2);
    mc_info_10_05.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_10_05.valid_bounds.max() = nvis::vec3(99., 99., 99.);
}




} // Garcia_vis_helper



#endif


































