#include <VTK/vtk_utils.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <list>
#include <image/nrrd_wrapper.hpp>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkWarpScalar.h>

struct lexi_less_with_tol {
    lexi_less_with_tol(double tol=1.0e-6) : _tol(tol) {}

    template<typename T, size_t N>
    bool operator()(const nvis::fixed_vector<T,N>& v0,
                    const nvis::fixed_vector<T,N>& v1) const {
        if (nvis::norm(v1-v0) < _tol) return false;
        else return nvis::lexicographical_order()(v0, v1);
    }

    double _tol;
};
typedef std::set<nvis::vec2, lexi_less_with_tol> set_type;
// notational convenience
typedef VTK_SMART(PolyData)      smart_poly;
typedef VTK_SMART(Actor)         smart_actor;
typedef VTK_SMART(Renderer)      smart_renderer;
typedef VTK_SMART(RenderWindow)  smart_window;

vtkActor* make_surface(const std::vector<nvis::vec2>& pts,
                       const std::vector<double>& vals,
                       double scale) {
    VTK_PTR(vtkPolyData, pd);
    pd = vtk_utils::make_points(pts);
    vtk_utils::add_scalars(pd, vals);
    vtk_utils::add_mesh2d(pd);
    VTK_PTR(vtkWarpScalar, warp);
    warp->SetInput(pd);
    warp->SetScaleFactor(scale);
    warp->Update();
    VTK_MAKE_ACTOR(an_actor, warp->GetOutput());
    an_actor->GetProperty()->SetColor(drand48(), drand48(), drand48());
    pd->Delete();
    warp->Delete();
    return an_actor;
}

int main(int argc, char* argv[]) {
    srand48(time(0));
    std::string filename = argv[1];
    double scale = 1;
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    std::vector<float> __array;
    xavier::to_vector(__array, nin);
    size_t nb_pts = __array.size()/5;

    set_type unique_src, unique_rec;
    std::list<smart_actor> surf_actors;

    nvis::bbox2 bounds;
    nvis::vec2 src, rec, last_src(0,0);
    std::vector<nvis::vec2> pts_on_surface;
    std::vector<double> ttime;
    for (size_t i=0 ; i<nb_pts ; ++i) {
        src[0] = __array[5*i  ];
        src[1] = __array[5*i+1];
        bounds.add(src);
        unique_src.insert(src);
        if (nvis::norm(src-last_src) > 1.0e-6) {
            if (i%10 == 1) { // visualize every 10th surface
                VTK_CREATE(vtkActor, sactor);
                sactor = make_surface(pts_on_surface, ttime, scale);
                surf_actors.push_back(sactor);
            }
            pts_on_surface.clear();
            ttime.clear();
            last_src = src;
        }
        rec[0] = __array[5*i+2];
        rec[1] = __array[5*i+3];
        bounds.add(rec);
        unique_rec.insert(rec);
        pts_on_surface.push_back(rec);
        ttime.push_back(__array[5*i+4]);
    }
    {
        VTK_CREATE(vtkActor, sactor);
        sactor = make_surface(pts_on_surface, ttime, scale);
        surf_actors.push_back(sactor);
    }
    nrrdNuke(nin);
    double sph_r = nvis::norm(bounds.size())/200.;

    std::vector<nvis::vec2> src_pts(unique_src.begin(), unique_src.end());
    std::vector<nvis::vec2> rec_pts(unique_rec.begin(), unique_rec.end());
    std::cout << "there are " << src_pts.size() << " sources and "
              << rec_pts.size() << " receivers\n";

    VTK_CREATE(vtkPolyData, src_pd);
    src_pd = vtk_utils::make_points(src_pts);
    vtk_utils::add_vertices(src_pd);
    // lift the sources
    VTK_CREATE(vtkTransform, trans);
    trans->Translate(0,0,10.*sph_r);
    VTK_CREATE(vtkTransformPolyDataFilter, tfilter);
    VTK_CONNECT(tfilter, src_pd);
    tfilter->SetTransform(trans);
    tfilter->Update();
    VTK_CREATE(vtkPolyData, src_pd_t);
    src_pd_t = tfilter->GetOutput();
    VTK_CREATE(vtkPolyData, rec_pd);
    rec_pd = vtk_utils::make_points(rec_pts);
    vtk_utils::add_vertices(rec_pd);
    VTK_CREATE(vtkPolyData, src_sph);
    src_sph = vtk_utils::make_spheres(src_pd_t, sph_r);
    VTK_CREATE(vtkPolyData, rec_sph);
    rec_sph = vtk_utils::make_spheres(rec_pd, sph_r);

    VTK_MAKE_ACTOR(src_actor, src_sph);
    VTK_MAKE_ACTOR(rec_actor, rec_sph);
    src_actor->GetProperty()->SetColor(1,1,0);
    rec_actor->GetProperty()->SetColor(0,0,1);

    VTK_CREATE(vtkRenderer, renderer);
    renderer->SetBackground(0,0,0);
    renderer->AddActor(src_actor);
    renderer->AddActor(rec_actor);
    for (std::list<smart_actor>::const_iterator it=surf_actors.begin() ;
         it!=surf_actors.end() ; ++it) {
        renderer->AddActor(*it);
    }
    VTK_CREATE(vtkRenderWindow, window);
    window->PointSmoothingOn();
    window->LineSmoothingOn();
    window->PolygonSmoothingOn();
    window->AddRenderer(renderer);
    window->SetSize(1200, 800);

    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    window->Render();
    interactor->Start();

    exit(0);
}