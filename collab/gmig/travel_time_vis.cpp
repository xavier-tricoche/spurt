#include <string>
#include <iostream>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkBMPReader.h>
#include <vtkCellArray.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourFilter.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGeometryFilter.h>
#include <vtkGlyph3D.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkMaskPoints.h>
#include <vtkPNGReader.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkTexture.h>
#include <vtkTIFFReader.h>
#include <vtkTubeFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkWarpScalar.h>

#include <VTK/vtk_utils.hpp>
#include <format/format.hpp>
#include "utils.hpp"
#include "travel_time.hpp"
#include "data_IO.hpp"

std::string me;
void usage(const std::string& msg)
{
    if (!msg.empty()) {
        std::cerr << "ERROR: " << msg << "\n\n";
    }
    std::cout
            << "USAGE: " << me << " [parameters] [options]\n"
            << '\n'
            << "DESCRIPTION: Visualize reconstructed travel time data\n"
            << '\n'
            << "PARAMETERS:\n"
            << " -i | --input <string>    Input smooth data reconstruction\n"
            << " -d | --data <string>     Input stations's data\n"
            << '\n'
            << "OPTIONS:\n"
            << " -h | --help              Print this information\n"
            << " -g | --gradient <string> Gradient file name            (default: \"\")\n"
            << " -x | --factor <float>    Scaling factor for 3D warping (default: 0)\n"
            << " -p | --path <string>     Path prepended to file names  (default: \"\")\n"
            << " -s | --show <float> x2   Show given position\n"
            << " -n | --number <int>      Number of isocontours to draw (default: 20)\n"
            << " -v | --verbose           Activate verbose mode         (default: false)\n"
            << "      --no-time           Don't show travel time isocontours\n"
            << "      --no-distance       Don't show distance to source isocontours\n"
            << "      --no-source         Don't show source location\n"
            << "      --no-stations       Don't show stations' locations\n";
    exit(1);
}

const double invalid_double = std::numeric_limits<double>::max();

// convert spherical coordinates to 3D position on unit sphere
inline nvis::vec3 to_sphere(const nvis::vec3& lola)
{
    // convert longitude and latitude to radians
    double phi = (lola[0]-90.)*M_PI/180.;
    double theta = lola[1]*M_PI/180.;
    return nvis::vec3(cos(phi)*cos(theta),
                      cos(phi)*sin(theta),
                      sin(phi));
}

inline double angular_distance(const nvis::vec3& x0, const nvis::vec3& x1)
{
    // project points on unit sphere to measure angle
    nvis::vec3 p0 = to_sphere(x0);
    nvis::vec3 p1 = to_sphere(x1);
    return acos(nvis::inner(p0, p1));
}

inline double geodesic_distance(const nvis::vec3& x0, const nvis::vec3& x1)
{
    const double Radius = 6371; // average earth radius in kilometers
    
    return Radius*angular_distance(x0, x1);
}

vtkColorTransferFunction*
blue_yellow_ctf(const vtkSmartPointer<vtkDataSet>& data)
{
    double range[2];
    data->GetPointData()->GetScalars()->GetRange(range);
    double min = range[0];
    double max = range[1];
    
    VTK_PTR(vtkColorTransferFunction, transfer_function);
    transfer_function->AddRGBPoint(min, 0, 0, 1);
    transfer_function->AddRGBPoint(0.75*min + 0.25*max, 0.5, 0.5, 1);
    transfer_function->AddRGBPoint(0.5*(min + max), 1, 1, 1);
    transfer_function->AddRGBPoint(0.25*min + 0.75*max, 1, 1, 0.5);
    transfer_function->AddRGBPoint(max, 1, 1, 0);
    return transfer_function;
}

vtkColorTransferFunction*
red_green_ctf(const vtkSmartPointer<vtkDataSet>& data)
{
    double range[2];
    data->GetPointData()->GetScalars()->GetRange(range);
    double min = range[0];
    double max = range[1];
    
    vtkColorTransferFunction* transfer_function = vtkColorTransferFunction::New();
    transfer_function->AddRGBPoint(min, 1, 0, 0);
    transfer_function->AddRGBPoint(0.75*min + 0.25*max, 1, 0.25, 0);
    transfer_function->AddRGBPoint(0.5*(min + max), 0.5, 0.75, 0);
    transfer_function->AddRGBPoint(0.25*min + 0.75*max, 0.25, 1, 0);
    transfer_function->AddRGBPoint(max, 0, 1, 0);
    return transfer_function;
}

int main(int argc, char* argv[])
{
    me = argv[0];
    std::string tt_name = "";
    std::string path = "";
    std::string station_name = "";
    std::string gradient_name = "";
    nvis::vec3 source;
    double scaling_factor = 0;
    bool relative_scaling = false;
    bool verbose = false;
    bool hide_time = false;
    bool hide_distance = false;
    bool hide_source = false;
    bool hide_stations = false;
    std::vector<nvis::vec2> show_pos;
    int nb_iso = 20;
    
    // parse command line parameters
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage("");
        } else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                usage("missing input file name");
            }
            tt_name = argv[++i];
        } else if (arg == "-d" || arg == "--data") {
            if (i == argc-1) {
                usage("missing data file name");
            }
            station_name = argv[++i];
        } else if (arg == "-g" || arg == "--gradient") {
            if (i == argc-1) {
                usage("missing gradient file name");
            }
            gradient_name = argv[++i];
        } else if (arg == "-n" || arg == "--number") {
           if (i == argc-1) {
               usage("missing number of isocontours");
           }
           nb_iso = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--show") {
            if (i >= argc-2) {
                usage("missing or incomplete position information");
            }
            show_pos.push_back(nvis::vec2());
            show_pos.back()[0] = atof(argv[++i]);
            show_pos.back()[1] = atof(argv[++i]);
        } else if (arg == "-x" || arg == "--factor") {
            if (i == argc-1) {
                usage("missing scaling factor");
            }
            std::string tmp_str = argv[++i];
            if (*tmp_str.rbegin() == '%') {
                scaling_factor = atof(tmp_str.substr(0, tmp_str.size()-1).c_str());
                relative_scaling = true;
            } else {
                scaling_factor = atof(tmp_str.c_str());
            }
        } else if (arg == "-p" || arg == "--path") {
            if (i == argc-1) {
                usage("missing path");
            }
            path = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "--no-time") {
            hide_time = true;
        } else if (arg == "--no-distance") {
            hide_distance = true;
        } else if (arg == "--no-source") {
            hide_source = true;
        } else if (arg == "--no-stations") {
            hide_stations = true;
        } else {
            usage("unrecognized argument");
        }
    }
    
    // check input
    
    if (source[0] == invalid_double) {
        usage("No source location provided");
    }
    if (tt_name.empty()) {
        usage("No input file name provided");
    }
    if (station_name.empty()) {
        usage("No data file name provided");
    }
    
    if (verbose && relative_scaling) {
        std::cout << "relative scaling factor for surface warp = " << scaling_factor << '\n';
    }
    
    // prepend file path, if available
    
    if (!path.empty()) {
        if (*path.rbegin() != '/') {
            path.append(1, '/');
        }
        tt_name = path + tt_name;
        if (!station_name.empty()) {
            station_name = path + station_name;
        }
    }
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //         Visualize travel time information as warped surface          //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    // import smoothly reconstructed travel time data
    vtkSmartPointer<vtkDataSet> tt_data;
    std::string tt_ext = xavier::get_extension(tt_name);
    if (tt_ext == "vtk") {
        VTK_CREATE(vtkDataSetReader, tt_reader);
        tt_reader->SetFileName(tt_name.c_str());
        tt_reader->Update();
        tt_data = tt_reader->GetOutput();
        tt_reader->GetOutput()->Register(tt_data);
        tt_reader->Delete();
    } else if (tt_ext == "nrrd" || tt_ext == "nhdr") {
        tt_data = vtk_utils::load_nrrd(tt_name);
    } else {
        usage("unrecognized file extension for " + tt_name);
    }
    tt_data->GetPointData()->GetScalars()->SetName("TravelTime");
    double _bounds[6];
    tt_data->ComputeBounds();
    tt_data->GetBounds(_bounds);
    nvis::vec3 _minp(_bounds[0], _bounds[2], _bounds[4]);
    nvis::vec3 _maxp(_bounds[1], _bounds[3], _bounds[5]);
    nvis::bbox3 bounds(_minp, _maxp);
    if (relative_scaling) {
        double _range[2];
        tt_data->GetPointData()->GetScalars()->GetRange(_range);
        double r = _range[1] - _range[0];
        scaling_factor *= nvis::norm(bounds.size())/r;
    }
    if (verbose) {
        std::cout << "spatial bounds of input dataset are " << bounds << '\n';
    }
    
    // assign texture coordinates to vertices for LIC visualization
    // of gradient information
    VTK_CREATE(vtkFloatArray, tcoords);
    int* dims = dynamic_cast<vtkImageData*>(tt_data.GetPointer())->GetDimensions();
    size_t npts = dims[0]*dims[1];
    tcoords->SetNumberOfComponents(3);
    tcoords->SetNumberOfTuples(npts);
    for (size_t n=0 ; n<npts ; ++n) {
        int i = n%dims[0];
        int j = n/dims[0];
        float tc[] = { (float)i/(float)(dims[0]-1), (float)j/(float)(dims[1]-1), 0 };
        tcoords->SetTuple(n, tc);
    }
    tt_data->GetPointData()->SetTCoords(tcoords);
    
    // create a color map for travel time
    VTK_CREATE(vtkColorTransferFunction, tt_tfx);
    tt_tfx = blue_yellow_ctf(tt_data);
    double* tt_range = tt_data->GetPointData()->GetScalars()->GetRange();
    double min_distance_to_source = tt_data->GetPointData()->GetScalars()->GetRange()[0];
    
    // warp travel time dataset to create heightfield representation
    VTK_CREATE(vtkWarpScalar, _warp);
    VTK_CONNECT(_warp, tt_data);
    _warp->SetScaleFactor(scaling_factor);
    
    // turn warped output into polydata for processing convenience
    VTK_CREATE(vtkGeometryFilter, warp);
    VTK_PLUG(warp, _warp);
    warp->Update();
    
    // render warped surface
    VTK_CREATE(vtkPolyDataMapper, tt_mapper);
    VTK_PLUG(tt_mapper, warp);
    tt_mapper->SetLookupTable(tt_tfx);
    VTK_CREATE(vtkActor, tt_actor);
    tt_actor->SetMapper(tt_mapper);
    tt_actor->GetProperty()->BackfaceCullingOff();
    tt_actor->GetProperty()->SetOpacity(1);
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //            Visualize station locations on warped surface             //
    //       (permits visual verification of interpolation property)        //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    // create sphere object to be used as a glyph
    VTK_CREATE(vtkSphereSource, sphere);
    sphere->SetThetaResolution(20);
    sphere->SetPhiResolution(20);
    sphere->SetRadius(0.005*nvis::norm(bounds.size()));
    
    // import station information from file
    xavier::gmig::traveltime::travel_time_data<double> station_tt_data;
    nvis::bbox2 b = xavier::gmig::traveltime::read_text(station_tt_data, 
                                                        station_name,
                                                        verbose);
    nvis::subv<0, 2, double, 3>(source) = station_tt_data.source;
    if (verbose) {
        std::cout << "source is located at " << source << '\n';
    }
    VTK_CREATE(vtkPolyData, stations);
    stations = vtk_utils::make_points(station_tt_data.receivers);
    stations->GetPoints()->GetData()->SetName("StationsLocations");
    vtk_utils::add_scalars(stations, station_tt_data.times);
    stations->GetPointData()->GetScalars()->SetName("StationsTravelTimes");
    vtk_utils::add_vertices(stations);
    
    // apply same warping to stations to lift them onto travel time surface
    VTK_CREATE(vtkWarpScalar, warp_s);
    VTK_CONNECT(warp_s, stations);
    warp_s->SetScaleFactor(scaling_factor);
    warp_s->Update();
    
    // apply glyph visualization to stations' locations
    VTK_CREATE(vtkGlyph3D, station_glyphs);
    VTK_PLUG(station_glyphs, warp_s);
    VTK_SOURCE_PLUG(station_glyphs, sphere);
    station_glyphs->ScalingOff();
    
    // render stations glyphs
    VTK_CREATE(vtkPolyDataMapper, station_mapper);
    station_mapper->ScalarVisibilityOff();
    VTK_PLUG(station_mapper, station_glyphs);
    VTK_CREATE(vtkActor, station_actor);
    station_actor->SetMapper(station_mapper);
    station_actor->GetProperty()->SetColor(1,0,1);
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //            Visualize user defined positions                          //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    VTK_CREATE(vtkActor, user_actor);
    if (show_pos.size()) {
        // import station information from file
        VTK_CREATE(vtkPolyData, user_pos);
        user_pos = vtk_utils::make_points(show_pos);
        user_pos->GetPoints()->GetData()->SetName("UserLocations");
        // assign travel time values to these positions
        std::vector<double> pos_time;
        for (int i=0 ; i<show_pos.size() ; ++i) {
            const nvis::vec2& x = show_pos[i];
            int nearest = tt_data->FindPoint(x[0], x[1], 0);
            if (nearest >= 0) {
                pos_time.push_back(tt_data->GetPointData()->GetScalars()->GetTuple1(nearest));
            }
        }
        vtk_utils::add_scalars(user_pos, pos_time);
        user_pos->GetPointData()->GetScalars()->SetName("UserLocationsTravelTimes");
        vtk_utils::add_vertices(user_pos);
        
        // apply same warping to stations to lift them onto travel time surface
        VTK_CREATE(vtkWarpScalar, warp_u);
        VTK_CONNECT(warp_u, user_pos);
        warp_u->SetScaleFactor(scaling_factor);
        warp_u->Update();
        
        // apply glyph visualization to stations' locations
        VTK_CREATE(vtkGlyph3D, user_glyphs);
        VTK_PLUG(user_glyphs, warp_u);
        VTK_SOURCE_PLUG(user_glyphs, sphere);
        user_glyphs->ScalingOff();
        
        // render stations glyphs
        VTK_CREATE(vtkPolyDataMapper, user_mapper);
        user_mapper->ScalarVisibilityOff();
        VTK_PLUG(user_mapper, user_glyphs);
        user_actor->SetMapper(user_mapper);
        user_actor->GetProperty()->SetColor(1,0,0);
    }
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //                     Visualize source location                        //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    // create a (trivial) dataset containing source as unique vertex
    std::vector<nvis::vec3> one_pt(1);
    one_pt[0] = source;
    one_pt[0][2] = scaling_factor*min_distance_to_source;
    std::vector<double> one_val(1);
    one_val[0] = 0.5;
    VTK_CREATE(vtkPolyData, source_pt);
    source_pt = vtk_utils::add_vertices(vtk_utils::make_points(one_pt));
    source_pt->GetPoints()->GetData()->SetName("SourceLocation");
    vtk_utils::add_scalars(source_pt, one_val);
    
    // apply glyph visualization to source dataset
    VTK_CREATE(vtkGlyph3D, source_sphere);
    VTK_CONNECT(source_sphere, source_pt);
    VTK_SOURCE_PLUG(source_sphere, sphere);
    source_sphere->ScalingOn();
    
    // render source
    VTK_CREATE(vtkPolyDataMapper, sphere_mapper);
    sphere_mapper->ScalarVisibilityOff();
    VTK_PLUG(sphere_mapper, source_sphere);
    VTK_CREATE(vtkActor, sphere_actor);
    sphere_actor->SetMapper(sphere_mapper);
    sphere_actor->GetProperty()->SetColor(1,1,0);
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //                Visualize isocontours of travel time                  //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    // compute isocontours
    VTK_CREATE(vtkContourFilter, tt_contour);
    VTK_PLUG(tt_contour, warp);
    double modified_tt_range[] = {tt_range[0] + 0.01*(tt_range[1]-tt_range[0]), tt_range[1]};
    tt_contour->GenerateValues(nb_iso, modified_tt_range);
    
    // wrap isocontours in tubes
    VTK_CREATE(vtkTubeFilter, tt_tube);
    VTK_PLUG(tt_tube, tt_contour);
    tt_tube->SetRadius(0.0025*nvis::norm(bounds.size()));
    tt_tube->SetNumberOfSides(20);
    
    // render travel time isocontours
    VTK_CREATE(vtkPolyDataMapper, tt_tube_mapper);
    VTK_PLUG(tt_tube_mapper, tt_tube);
    tt_tube_mapper->ScalarVisibilityOff();
    VTK_CREATE(vtkActor, tt_tube_actor);
    tt_tube_actor->SetMapper(tt_tube_mapper);
    tt_tube_actor->GetProperty()->SetColor(1,0,0);
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //             Visualize isocontours of distance to source              //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    // compute distances to source at each vertex
    size_t nb_points = tt_data->GetNumberOfPoints();
    std::vector<double> __distances(nb_points);
    for (size_t i=0 ; i<nb_points ; ++i) {
        nvis::vec3 x;
        tt_data->GetPoint(i, x.begin());
        double d = geodesic_distance(x, source);
        __distances[i] = d;
    }
    
    // assign distance information to warped surface
    VTK_CREATE(vtkPolyData, distance);
    distance->CopyStructure(warp->GetOutput());
    vtk_utils::add_scalars(distance, __distances);
    distance->GetPointData()->GetScalars()->SetName("DistanceToSource");
    double dist_range[2];
    distance->GetPointData()->GetScalars()->GetRange(dist_range);
    
    // create a color map for distance data
    VTK_CREATE(vtkColorTransferFunction, dist_tfx);
    dist_tfx = red_green_ctf(distance);
    
    // render distance data on warped surface
    VTK_CREATE(vtkPolyDataMapper, dist_mapper);
    VTK_CONNECT(dist_mapper, distance);
    dist_mapper->SetLookupTable(dist_tfx);
    VTK_CREATE(vtkActor, dist_actor);
    dist_actor->SetMapper(dist_mapper);
    
    // compute isocontours in distance field
    VTK_CREATE(vtkContourFilter, dist_contour);
    VTK_CONNECT(dist_contour, distance);
    dist_contour->GenerateValues(nb_iso, dist_range);
    double* values = dist_contour->GetValues();
    dist_contour->Update();
    
    // wrap distance isocontours in tubes
    VTK_CREATE(vtkTubeFilter, dist_tube);
    VTK_PLUG(dist_tube, dist_contour);
    dist_tube->SetRadius(0.0025*nvis::norm(bounds.size()));
    dist_tube->SetNumberOfSides(20);
    
    // render distance iscontours
    VTK_CREATE(vtkPolyDataMapper, dist_tube_mapper);
    VTK_PLUG(dist_tube_mapper, dist_tube);
    dist_tube_mapper->ScalarVisibilityOff();
    VTK_CREATE(vtkActor, dist_tube_actor);
    dist_tube_actor->SetMapper(dist_tube_mapper);
    dist_tube_actor->GetProperty()->SetColor(0,1,0);
    
    // create rendering engine and render window
    VTK_CREATE(vtkRenderer, renderer);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(800, 800);
    
    // -------------------------------------------------------------------- //
    //                                                                      //
    //     Visualize gradient information as LIC texture, if available      //
    //                                                                      //
    // -------------------------------------------------------------------- //
    
    VTK_CREATE(vtkActor, lic_actor);
    if (!gradient_name.empty()) {
        vtkImageData* lic_image;
        // check if the user already provided us with an image,
        // in which case we will assume that it corresponds to a
        // precomputed LIC texture.
        std::string ext = xavier::get_extension(gradient_name);
        if (ext == "png" || ext == "jpeg" || ext == "jpg" || ext == "jpg2" ||
                ext == "bmp" || ext == "jpeg2" || ext == "tiff" || ext == "tif") {
            lic_image = vtk_utils::load_image(gradient_name);
        } else {
            // import gradient vector data from file
            vtkImageData* rhs;
            if (ext == "nrrd" || ext == "nhdr") {
                rhs = vtk_utils::load_nrrd(gradient_name);
            } else {
                VTK_PTR(vtkStructuredPointsReader, reader);
                reader->SetFileName(gradient_name.c_str());
                reader->Update();
                rhs = static_cast<vtkImageData*>(reader->GetOutput());
                rhs->Register(rhs);
                reader->Delete();
            }
            // compute LIC texture
            lic_image = vtk_utils::do_lic(rhs, 100, 0.5, 4);
            rhs->Delete();
        }
        // assign LIC image to texture
        vtkTexture* lic_texture = vtk_utils::image_to_texture(lic_image);
        lic_image->Delete();
        // assign texture to surface actor
        tt_actor->SetTexture(lic_texture);
        lic_texture->Delete();
        // since we have a texture to display, turn off travel time color map
        tt_actor->GetMapper()->ScalarVisibilityOff();
    }
    
    // add actors to rendering engine
    renderer->AddActor(tt_actor);             // * travel time surface
    if (!hide_time) {
        renderer->AddActor(tt_tube_actor);    // * travel time isoncontours
    }
    if (!hide_distance) {
        renderer->AddActor(dist_tube_actor);  // * distance to source isocontours
    }
    if (!hide_source) {
        renderer->AddActor(sphere_actor);     // * source location
    }
    if (!hide_stations) {
        renderer->AddActor(station_actor);    // * stations locations
    }
    if (show_pos.size()) {
        renderer->AddActor(user_actor);       // * user selected locations
    }
    
    // configure rendering window
    renderer->SetBackground(0,0,0);
    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    renderer->LightFollowCameraOn();
    
    // enter interaction loop
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    window->Render();
    interactor->Initialize();
    interactor->Start();
    
    return 0;
}