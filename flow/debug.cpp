#include <fstream>
#include <iostream>
#include <string>
#include <map>

#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>

#include <math/fixed_vector.hpp>

#include <VTK/vtk_utils.hpp>

#include <vtkTransform.h>
#include <vtkTransformFilter.h>

const std::string mask_name = "/home/xmt/visdata/Flows/CFD/NCOM/nrrd/land_mask.nrrd";

const float LAND = -1000;
const float SEA = -500;
float radius = 0.001;
bool notubes = false;
double min_dist = 1.0e-6;
float time_scaling = 0;
std::vector<int> only_indices;

std::string name_in, name_out;
std::vector<std::string> input_files;
typedef nvis::fixed_vector<double, 4> vec4;
typedef nvis::fixed_vector<double, 6> vec6;
vec4 bnds;

template< typename T >
inline T nrrd_value(const Nrrd* nin, size_t n) {
    return spurt::nrrd_data_wrapper<T>(nin)[n];
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Visualize results of particle integration");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_sequence("input", input_files, "Input trajectory filenames", required_group);
        parser.add_value("scale", time_scaling, time_scaling, "Scaling coefficient for time axis", optional_group);
        parser.add_value("radius", radius, radius, "tube radius", optional_group);
        parser.add_value("output", name_out, "Output file basename", required_group);
        parser.add_flag("notubes", notubes, "Do not use tubes for display", optional_group);
        parser.add_sequence("index", only_indices, "Display only selected trajectory indices", optional_group);
        parser.add_tuple<4>("bounds", bnds, bnds, "Sampling bounds (min longitude and latitude followed by max's)", required_group);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

VTK_SMART(Actor) make_tubes(vtkPolyData* input, const nvis::vec3& col) {
    VTK_CREATE(vtkPolyDataMapper, mapper);
    VTK_CREATE(vtkActor, actor);
    if (!notubes) {
        VTK_CREATE(vtkTubeFilter, tubes);
        tubes->SetInputData(input);
        tubes->SetRadius(radius);
        tubes->SetNumberOfSides(6);
        mapper->SetInputConnection(tubes->GetOutputPort());
    }
    else {
        mapper->SetInputData(input);
    }
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(col[0], col[1], col[2]);
    return actor;
}

int main(int argc, char* argv[]) {
    initialize(argc, const_cast<const char**>(argv));
    
    VTK_CREATE(vtkRenderer, renderer);
    
    Nrrd* mask = spurt::readNrrd(mask_name);
    nvis::bbox2 region;
    region.min() = nvis::vec2(bnds[0], bnds[1]);
    region.max() = nvis::vec2(bnds[2], bnds[3]);
    nvis::bbox2 bounds;
    nvis::ivec2 res(mask->axis[0].size, mask->axis[1].size);
    nvis::vec2 spc(mask->axis[0].spacing, mask->axis[1].spacing);
    nvis::vec2 mins(mask->axis[0].min, mask->axis[1].min);
    nvis::vec2 maxs(mins[0] + (res[0]-1)*spc[0], mins[1] + (res[1]-1)*spc[1]);
    
    std::vector<nvis::vec2> region_pts(4);
    region_pts[0] = region.min();
    region_pts[1] = nvis::vec2(region.max()[0], region.min()[1]);
    region_pts[2] = region.max();
    region_pts[3] = nvis::vec2(region.min()[0], region.max()[1]);
    vtkSmartPointer<vtkPolyData> region_pd = vtk_utils::make_points(region_pts);
        
    std::vector<int> region_tris(6);
    region_tris[0] = 0;
    region_tris[1] = 1;
    region_tris[2] = 2;
    region_tris[3] = 0;
    region_tris[4] = 2;
    region_tris[5] = 3;
    vtk_utils::add_mesh2d(region_pd, region_tris);
    
    VTK_CREATE(vtkTransformFilter, tfilter);
    VTK_CREATE(vtkTransform, transform);
    transform->Translate(0,0,-0.5);
    tfilter->SetTransform(transform.GetPointer());
    tfilter->SetInputData(region_pd);
    tfilter->Update();
    
    VTK_MAKE_ACTOR(region_actor, tfilter->GetOutput());
    region_actor->GetProperty()->SetColor(1,1,0);
    region_actor->GetProperty()->SetOpacity(0.5);
    // renderer->AddActor(region_actor);
    
    VTK_CREATE(vtkStructuredPoints, domain);
    domain->SetDimensions(res[0], res[1], 1);
    domain->SetOrigin(mins[0], mins[1], -1);
    domain->SetSpacing(spc[0], spc[1], 1);
    std::vector<float> values;
    spurt::to_vector(values, mask);
    vtk_utils::add_scalars(domain, values);
    VTK_CREATE(vtkDataSetMapper, domain_mapper);
    domain_mapper->SetInputData(domain);
    domain_mapper->ScalarVisibilityOn();
    VTK_CREATE(vtkColorTransferFunction, cmap1);
    cmap1->AddRGBPoint(0, 0, 0.75, 0);
    cmap1->AddRGBPoint(1, 0, 0.2, 1);
    domain_mapper->SetLookupTable(cmap1);
    VTK_CREATE(vtkActor, domain_actor);
    domain_actor->SetMapper(domain_mapper);
    renderer->AddActor(domain_actor);
    
    std::vector<nvis::vec2> seeds;
    std::vector<nvis::vec2> samples;
    
    std::map<int, int> id_to_id;
    std::vector< std::vector<nvis::vec3> > pathlines;
    std::vector< std::vector<nvis::vec2> > streamlines;
    std::vector< std::vector< int > > ids;
    std::vector< std::vector< float > > lavds;
    
    std::set<int> indices;
    if (!only_indices.empty()) {
        std::copy(only_indices.begin(), only_indices.end(), std::inserter(indices, indices.begin()));
    }
    
    for (size_t i=0; i<input_files.size(); ++i) {
        name_in = input_files[i];
        spurt::ProgressDisplay progress(false);
        Nrrd* traj_nrrd = spurt::readNrrd(name_in);
        float* data = (float*)traj_nrrd->data;
        size_t n_pts = traj_nrrd->axis[1].size;
        progress.start(n_pts, "Processing filename #" + std::to_string(i));
        progress.active = true;
        int id;
        for (size_t n=0; n<n_pts ; ++n) {
            id = static_cast<int>(data[5*n+4]);
            if (!indices.empty() && indices.find(id)==indices.end()) continue;
            float x = data[5*n];
            float y = data[5*n+1];
            float t = data[5*n+2];
            if (time_scaling > 0) t*=time_scaling;
            float lavd = data[5*n+3];
            auto it = id_to_id.find(id);
            if (it == id_to_id.end()) {
                pathlines.push_back(std::vector<nvis::vec3>());
                ids.push_back(std::vector< int >());
                lavds.push_back(std::vector< float >());
                id_to_id[id] = pathlines.size()-1;
                pathlines.back().push_back(nvis::vec3(x,y,t));
                ids.back().push_back(id);
                lavds.back().push_back(lavd);
            }
            else {
                pathlines[it->second].push_back(nvis::vec3(x,y,t));
                ids[it->second].push_back(id);
                lavds[it->second].push_back(lavd);
            }
            progress.update(n);
        }
        progress.end();
        std::cout << "Compute time: cpu: " << progress.cpu_time() << " ms.) |"
            << " wall: " << progress.wall_time() << " ms.)\n";
    
        std::cout << "max line id = " << id << '\n';
        nrrdNuke(traj_nrrd);
    }
    
    struct time_sort {
        bool operator()(const nvis::vec3& a, const nvis::vec3& b) {
            return a[2] > b[2];
        }
    };
    
    size_t mean_length=0;
    for (size_t n=0; n<pathlines.size(); ++n) {
        if (pathlines[n].empty()) continue;
        std::sort(pathlines[n].begin(), pathlines[n].end(), time_sort());
        if (time_scaling == 0) {
            streamlines.push_back(std::vector<nvis::vec2>());
            for (size_t i=0; i<pathlines[n].size(); ++i) {
                streamlines.back().
                    push_back(nvis::subv<0, 2, double, 3>(pathlines[n][i]));
            }
        }
        mean_length += pathlines[n].size();
    }
    std::cout << "Average pathline length=" << mean_length/(float)pathlines.size() << '\n';
    
    std::vector<nvis::ivec2> removed;
    
    vtkSmartPointer<vtkPolyData> sl_data;
    if (!time_scaling)
        sl_data = vtk_utils::make_polylines(streamlines, removed, min_dist);
    else 
       sl_data = vtk_utils::make_polylines(pathlines, removed, min_dist);
    std::cout << removed.size() << " removed vertices\n";
    
    std::set<nvis::ivec2, nvis::lexicographical_order> missing;
    std::copy(removed.begin(), removed.end(), std::inserter(missing, missing.begin()));

    std::vector<int> clean_ids;
    std::vector<float> clean_lavds;
    for (size_t n=0; n<ids.size(); ++n) {
        for (size_t i=0; i<ids[n].size(); ++i) {
            if (missing.find(nvis::ivec2(n,i)) == missing.end()) {
                clean_ids.push_back(ids[n][i]);
                clean_lavds.push_back(lavds[n][i]);
            }
        }
    }
    std::cout << sl_data->GetNumberOfPoints() << " points in sl_data\n";
    std::cout << clean_ids.size() << " values in ids\n";
    std::cout << clean_lavds.size() << " values in lavds\n";

    
    int max_id = *std::max_element(clean_ids.begin(), clean_ids.end());
    float min_lavd, max_lavd;
    min_lavd = *std::min_element(clean_lavds.begin(), clean_lavds.end());
    max_lavd = *std::max_element(clean_lavds.begin(), clean_lavds.end());
    
    std::cout << "min lavd = " << min_lavd << ", max lavd = " << max_lavd << '\n';
    
    vtk_utils::add_scalars(sl_data, clean_lavds);
    
    std::sort(clean_lavds.begin(), clean_lavds.end());
    size_t stride = clean_lavds.size()/10;
    const nvis::fvec3 yellow(1,1,0);
    const nvis::fvec3 blue(0,0,1);
    
    VTK_CREATE(vtkColorTransferFunction, ctf);
    for (int i=0; i<clean_lavds.size(); i+=stride) {
        float u = i/(float)clean_lavds.size();
        nvis::fvec3 c = (1.-u)*blue + u*yellow;
        ctf->AddRGBPoint(clean_lavds[i], c[0], c[1], c[2]);
    }
    VTK_SMART(Actor) sl_actor = make_tubes(sl_data, nvis::vec3(0, 0.5, 1));
    sl_actor->GetMapper()->ScalarVisibilityOn();
    sl_actor->GetMapper()->SetLookupTable(ctf);
    renderer->AddActor(sl_actor);
    
    vtk_utils::colorbar_param cbparam;
    cbparam.height=300;
    cbparam.pos[0] = 0.85;
    cbparam.title = "LAVD";
    
    VTK_SMART(ScalarBarActor) scalarbar = vtk_utils::colorbar(ctf, cbparam);
    scalarbar->SetNumberOfLabels(10);
    
    /*
    std::vector< std::vector<nvis::vec2> > failed_paths;
    std::vector<float> failed_paths_ids;
    progress.start(failed.size(), "Processing failed paths");
    progress.active = true;
    for (size_t n=0; n<failed.size() && n<1; ++n) {
        if (drand48()>0.1) continue;
        failed_paths.push_back(std::vector<nvis::vec2>());
        std::vector<nvis::vec2>& sl = failed_paths.back();
        Nrrd* nin = spurt::readNrrd(failed[n]);
        size_t npts = nin->axis[1].size;
        std::cout << "\n" << npts << " points in failed path\n";
        for (int k=0; k<npts; ++k) {
            float x = nrrd_value<float>(nin, 4*k);
            float y = nrrd_value<float>(nin, 4*k+1);
            sl.push_back(nvis::vec2(x,y));
            failed_paths_ids.push_back(cos(n/10.));
        }
        progress.update(n);
    }
    progress.end();
    std::cout << "Compute time: cpu: " << progress.cpu_time() << " ms.) |"
        << " wall: " << progress.wall_time() << " ms.)\n";
    
    vtkSmartPointer<vtkPolyData> failed_data = vtk_utils::make_polylines(failed_paths);
    vtk_utils::add_scalars(failed_data, failed_paths_ids);
    VTK_CREATE(vtkColorTransferFunction, cmap2);
    cmap2->AddRGBPoint(-1, 1, 1, 1);
    cmap2->AddRGBPoint(0, 1, 0.5, 0.5);
    cmap2->AddRGBPoint(1, 1, 0, 0);    
    VTK_SMART(Actor) failed_actor = make_tubes(failed_data, nvis::vec3(1,0,0));
    renderer->AddActor(failed_actor);
    */
    
    renderer->AddActor2D(scalarbar);
    renderer->SetBackground(0, 0, 0);
    renderer->ResetCamera();
    
    nrrdNuke(mask);
    
    VTK_CREATE(vtkRenderWindow, window);
    window->SetSize(1024, 768);
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    window->AddRenderer(renderer);
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    window->Render();
    interactor->Start();
    
    return 0;
}
