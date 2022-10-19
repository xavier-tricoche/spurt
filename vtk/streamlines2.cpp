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
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkScalarBarActor.h"
#include "vtkColorTransferFunction.h"
#include "vtkTextProperty.h"

#include <string>
#include <math/fixed_vector.hpp>
#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <teem/hest_helper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include <vis/streamline.hpp>
#include <boost/shared_ptr.hpp>

inline void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

char*    file;
char*    scalar;
int     nblines;
float   length;
float   center[3];
float   radius[3];
float   eps;
float   min_step;
int     discretization;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &file,                NULL,       "input vector field (NRRD or VTK)");
    hestOptAdd(&hopt, "s",      "scalar",           airTypeString,  0, 1, &scalar,              NULL,       "scalar field used for color coding (NRRD)");
    hestOptAdd(&hopt, "n",      "nb lines",         airTypeInt,     0, 1, &nblines,             "100",      "number of lines");
    hestOptAdd(&hopt, "l",      "length",           airTypeFloat,   1, 1, &length,              NULL,       "integration length");
    hestOptAdd(&hopt, "c",      "center",           airTypeFloat,   3, 3, &center,              NULL,       "seed center");
    hestOptAdd(&hopt, "r",      "radius",           airTypeFloat,   3, 3, &radius,              NULL,       "seed radius");
    hestOptAdd(&hopt, "e",      "eps",              airTypeFloat,   0, 1, &eps,                 "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "min",    "min step",         airTypeFloat,   0, 1, &min_step,            "1.0e-8",   "step size underflow threshold");
    hestOptAdd(&hopt, "d",      "discretization",   airTypeInt,     0, 1, &discretization,      "20",       "number of lines");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Compute and visualize streamlines in 3D vector field",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}



template<typename T>
struct field {

    typedef T       value_type;
    
    template<typename T1>
    field(const std::vector<T1>& data,
          const nvis::ivec3& size,
          const nvis::bbox3& bounds)
        : __data(data.size()), __size(size), __bounds(bounds) {
        for (int i = 0 ; i < sz ; ++i) {
            __data[i] = data[i];
        }
        __step = __bounds.size() / nvis::vec3(__size - nvis::ivec3(1, 1, 1));
    }
    
    template<typename T1>
    field(const Nrrd* nin) {
        for (int i = 0 ; i < 3 ; ++i) {
            __size[i] = nin->axis[nin->dim-3+i].size;
        }
        size_t size = __size[0] * __size[1] * __size[2];
        __bounds = xavier::nrrd_utils::get_bounds(nin);
        __step = __bounds.size() / nvis::vec3(__size - nvis::ivec3(1, 1, 1));
        __data = std::vector<value_type>((T1*)nin->data, (T1*)nin->data + size);
    }
    
    vector_field(const vtkSmartPointer<vtkDataSet>& dataset)
        : __data(dataset->GetNumberOfPoints()) {
        if (__vtk_data->GetDataDimension() != 3) {
            throw;
        }
        
        dataset->GetDimensions(__size.begin());
        double bounds[6];
        dataset->GetBounds(bounds);
        __bounds.min() = nvis::vec3(bounds[0], bounds[2], bounds[4]);
        __bounds.max() = nvis::vec3(bounds[1], bounds[3], bounds[5]);
        __step = __bounds.size() / nvis::vec3(__size - nvis::ivec3(1, 1, 1));
        
        double* val;
        for (int i = 0 ; i < __data.size() ; ++i) {
            val = dataset->GetPointData()->GetVectors()->GetTuple(i);
            __data[i] = value_type(val[0], val[1], val[2]);
        }
    }
    
    ~vector_field() {}
    
    bool g2l(nvis::vec3& l, const nvis::vec3& g) const {
        if (!__bounds.inside(g)) {
            return false;
        }
        l = (g - __bounds.min()) / __step;
        return true;
    }
    
    int index(const nvis::ivec3& id) const {
        return id[0] + __size[0]*(id[1] + __size[1]*id[2]);
    }
    
    bool operator()(double, const nvis::vec3& x, value_type& v) const {
        nvis::vec3 y;
        if (!g2l(y, x)) {
            return false;
        }
        nvis::ivec3 id(floor(y[0]), floor(y[1]), floor(y[2]));
        nvis::vec3 u = y - nvis::vec3(id);
        double w[] = {(1. - u[0])* (1. - u[1])* (1. - u[2]),
                      u[0]* (1. - u[1])* (1. - u[2]),
                      u[0]* u[1]* (1. - u[2]),
                      (1. - u[0])* u[1]* (1. - u[2]),
                      (1. - u[0])* (1. - u[1])* u[2],
                      u[0]* (1. - u[1])* u[2],
                      u[0]* u[1]* u[2],
                      (1. - u[0])* u[1]* u[2]
                     };
        v = w[0] * __data[index(id)] +
            w[1] * __data[index(id+nvis::ivec3(1, 0, 0))] +
            w[2] * __data[index(id+nvis::ivec3(1, 1, 0))] +
            w[3] * __data[index(id+nvis::ivec3(0, 1, 0))] +
            w[4] * __data[index(id+nvis::ivec3(0, 0, 1))] +
            w[5] * __data[index(id+nvis::ivec3(1, 0, 1))] +
            w[6] * __data[index(id+nvis::ivec3(1, 1, 1))] +
            w[7] * __data[index(id+nvis::ivec3(0, 1, 1))];
            
        return true;
    }
    
    std::vector<value_type>     __data;
    nvis::ivec3                 __size;
    nvis::bbox3                 __bounds;
    nvis::vec3                  __step;
};

inline vtkSmartPointer<vtkPolyData> to_polydata(const std::list<nvis::vec3>& line)
{
    vtkSmartPointer<vtkDoubleArray> coords = vtkSmartPointer<vtkDoubleArray>::New();
    coords->SetNumberOfComponents(3);
    for (std::list<nvis::vec3>::const_iterator it = line.begin() ;
            it != line.end() ; ++it) {
        coords->InsertNextTuple(it->begin());
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coords);
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType n = 2;
    for (int i = 0 ; i < line.size() - 1 ; ++i) {
        vtkIdType ids[] = { i, i + 1 };
        lines->InsertNextCell(n, ids);
    }
    
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->SetLines(lines);
    return pd;
}

struct stop {
    enum {
        LEFT_DOMAIN,
        CRITICAL_POINT,
    } state;
    
    stop(double eps = 0.) : __eps(eps) {}
    
    bool operator()(const nvis::streamline::int_step& step) {
        if (step.y1()[4] != 0) {
            state = LEFT_DOMAIN;
            return true;
        } else if (fabs(step.y0()[3] - step.y1()[3]) < __eps) {
            state = CRITICAL_POINT;
            return true;
        } else {
            return false;
        }
    }
    
    double __eps;
};

std::ostream& operator<<(std::ostream& os, const stop& __stop)
{
    if (__stop.state == stop::LEFT_DOMAIN) {
        os << "integration left domain";
    } else if (__stop.state == stop::CRITICAL_POINT) {
        os << "integration reached critical point";
    } else {
        os << "integration successful\n";
    }
    return os;
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    typedef nvis::fvec3         value_type;
    typedef field<value_type>   field_type;
    
    boost::shared_ptr<field_type> vf;
    
    // determine extension
    std::string name(file);
    if (name.substr(name.find_last_of(".") + 1) == "nrrd") {
        std::cerr << "NRRD file extension recognized.\n";
        Nrrd* nin = xavier::nrrd_utils::readNrrd(file);
        nvis::ivec3 size;
        if (nin->dim != 4) {
            throw;
        }
        for (int i = 0 ; i < 3 ; ++i) {
            size[i] = nin->axis[i+1].size;
        }
        int sz = size[0] * size[1] * size[2];
        value_type* vdata = new value_type[sz];
        std::vector<float> sdata;
        xavier::to_vector(sdata, nin);
        for (int i = 0 ; i < sz ; ++i) {
            data[i][0] = sdata[3*i  ];
            data[i][1] = sdata[3*i+1];
            data[i][2] = sdata[3*i+2];
        }
        vf = boost::shared_ptr<field_type>(new field_type(data, size, xavier::nrrd_utils::get_bounds(nin)));
        nrrdNuke(nin);
    } else if (name.substr(name.find_last_of(".") + 1) == "vtk") {
        std::cerr << "VTK file extension recognized.\n";
        vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
        reader->SetFileName(file);
        reader->Update();
        vtkStructuredPoints* data = reader->GetOutput();
        vf = boost::shared_ptr<vector_field>(new vector_field(data));
    } else {
        std::cerr << "file extension in " << name << " was not recognized.\n";
        return 1;
    }
    
    bool use_norm = (scalar == NULL);
    if (!use_norm) {
    
    }
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0, 0, 0);
    
    std::cerr << "Camera set\n";
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    // ren_win->PointSmoothingOn();
    // ren_win->LineSmoothingOn();
    // ren_win->PolygonSmoothingOn();
    ren_win->SetAlphaBitPlanes(1);
    // ren_win->SetMultiSamples(10);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    
    nvis::vec3 min(center[0] - radius[0],
                   center[1] - radius[1],
                   center[2] - radius[2]);
    nvis::vec3 r(2.*radius[0], 2.*radius[1], 2.*radius[2]);
    
    srand48(0);
    for (int i = 0 ; i < nblines ; ++i) {
        std::list<nvis::vec3> line;
        
        nvis::vec3 seed = min + nvis::vec3(drand48(), drand48(), drand48()) * r;
        
        std::cerr << "seeding streamline #" << i << "/" << nblines << " at " << seed << '\n';
        
        nvis::streamline sl(seed);
        sl.record = true;
        sl.reltol = sl.abstol = eps;
        sl.stepsz = 1.;
        
        stop stop_criterion(min_step);
        int fwd = sl.advance(*vf, length, stop_criterion);
        if (fwd == nvis::streamline::CUSTOM) {
            std::cerr << "forward integration stopped because " << stop_criterion << '\n';
        }
        int bwd = sl.advance(*vf, -length, stop_criterion);
        if (bwd == nvis::streamline::CUSTOM) {
            std::cerr << "backward integration stopped because " << stop_criterion << '\n';
        }
        
        double dt = (sl.t_max() - sl.t_min()) / (float)discretization;
        for (double t = sl.t_min() ; t <= sl.t_max() ; t += dt) {
            line.push_back(sl(t));
            // std::cerr << "added position sl(" << t << ") = " << line.back() << '\n';
        }
        vtkSmartPointer<vtkPolyData> pd = to_polydata(line);
        
        vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
        tubes->SetInput(pd);
        tubes->SetRadius(1);
        tubes->SetNumberOfSides(6);
        
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(tubes->GetOutputPort());
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1, 1, 1);
        
        ren->AddActor(actor);
    }
    
    std::cerr << "done.\n";
    
    ren->ResetCamera();
    
    ren_win->Render();
    iren->Initialize();
    iren->Start();
    
    return 0;
}































































