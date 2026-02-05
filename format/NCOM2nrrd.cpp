#include <numeric>
#include <array>

#include <format/NCOMreader.hpp>
#include <format/filename.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/strings.hpp>
#include <misc/progress.hpp>

#include <vtk/vtk_data_helper.hpp>
#include <third_party/spline/spline.h>

std::string name_in, name_out, var_names;
bool verbose=false;
bool save_mesh=false;
bool save_uv=true;
double ref_lat=24.42, ref_lon=270.9;
std::string kernel_name = "Scharr";

constexpr double D2R = M_PI/180.;
constexpr double R2D = 180./M_PI;
constexpr double Earth_radius = 6371000; // in meters
const double _invalid_ = -30000;

inline double deg2rad(double d) {
    return d*D2R;
}

inline double rad2deg(double r) {
    return r*R2D;
}

/*
    polar model of an ellipse:
        x, y = a*cos(t), b*sin(t)
    
    t *is not* the latitude in this case. Instead, it satisfies the following
    equation:
        lat(t) = atan(y/x) = atan(asin(t)/bcos(t)) = atan(a/b*tan(t))

    The radius of parallel at t is a*cos(t). It is the radius of parallel at
    latitude atan(a/b*tan(t))
    The radius of the ellipse at t is sqrt((a*cos(t))**2 + (b*sin(t))**2)

    We create a table of t <-> lat pairs that we interpolate to find the t value 
    for a given latitude such that we can plug it in the expressions above.
*/

struct oblate {
    static constexpr double a=6378137;
    static constexpr double b=6356752;

    tk::spline m_spline;

    oblate(int res=100) : m_spline() {
        std::vector<double> ts(res);
        std::iota(ts.begin(), ts.end(), 0);
        double step = M_PI/2./(res-1.);
        std::for_each(ts.begin(), ts.end(), [&](double& t){ t*=step; });
        std::vector<double> thetas(ts.begin(), ts.end());
        std::for_each(thetas.begin(), thetas.end(), [&](double& t){ t=std::atan(a/b*std::tan(t)); });
        m_spline.set_points(thetas, ts); // interpolates theta -> t mappings
    }

    double radius(double lat) {
        double t = m_spline(lat); // interpolate to get t for given lat
        return std::sqrt(a*a*std::cos(t)*std::cos(t) + b*b*std::sin(t)*std::sin(t));
    }

    double parallel_radius(double lat) {
        double t = m_spline(lat); // interpolate to get t for given lat
        return a*std::cos(t);
    }
};

using namespace spurt;

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Convert NCOM file to NRRD format");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required_group);
        parser.add_value("output", name_out, "Output filename", optional_group);
        parser.add_flag("mesh", save_mesh, "Export mesh information", optional_group);
        parser.add_flag("uv", save_uv, "Export velocity in m/s", optional_group);
        parser.add_value("vars", var_names, "Variable names", optional_group);
        parser.add_value("kernel", kernel_name, "Kernel name", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);
        
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

void check_delta(const std::vector<double>& array) {
    std::vector<double> deltas;
    for (size_t i=0; i<array.size()-1; ++i) {
        deltas.push_back(array[i+1]-array[i]);
    }
/*    
 *  double mean=std::accumulate(deltas.begin(), deltas.end(), 0)/deltas.size();
    double min=*std::min_element(deltas.begin(), deltas.end());
    double max=*std::max_element(deltas.begin(), deltas.end());
    
    std::cout << std::setprecision(12) << "mean step=" << mean << ", min step=" << min << ", max step=" << max << '\n'; 
    */
}

int main(int argc, char* argv[]) {
    initialize(argc, (const char**)argv);
    
    std::vector<double> lat_deg, lon_deg, lat_rad, lon_rad;
    std::vector<vec2> vel_spatial, vel_angular;
    size_t nlat, nlon;
    double t;
    if (name_out.empty()) {
        name_out = spurt::filename::remove_extension(name_in);
    }
    else {
        name_out = spurt::filename::remove_extension(name_out);
    }
    
    std::vector<std::string> names;
    std::vector< std::vector<double > > variables;
    if (!var_names.empty()) {
        spurt::tokenize(names, var_names);
        spurt::NCOMreader::load_dataset(name_in, lat_deg, lon_deg, vel_spatial, names, variables, t);
    }
    else {
        spurt::NCOMreader::load_dataset(name_in, lat_deg, lon_deg, vel_spatial, t);
    }
    nlat=lat_deg.size();
    nlon=lon_deg.size();
    if (verbose) {
        std::cout << "latitude range (in degrees): " << lat_deg[0] << " - " << lat_deg.back() << '\n';
        std::cout << "longitude range (in degrees): " << lon_deg[0] << " - " << lon_deg.back() << '\n';
    }
    
    std::vector<double> vel_norm;
    size_t nvalid = 0;
    for (size_t i=0; i<vel_spatial.size(); ++i) {
        if (vel_spatial[i][0] == _invalid_ || vel_spatial[i][1] == _invalid_) {
            vel_spatial[i] = vec2(0,0);
        }
        else {
            vel_spatial[i] *= 0.001;
            vel_norm.push_back(norm(vel_spatial[i]));
        }
    }
    
    double max_vel, min_vel, mean_vel=0;
    max_vel = *std::max_element(vel_norm.begin(), vel_norm.end());
    min_vel = *std::min_element(vel_norm.begin(), vel_norm.end());
    std::for_each(vel_norm.begin(), vel_norm.end(), [&](double v) 
    {
        mean_vel += v;
    });
    mean_vel /= vel_norm.size();
    
    if (verbose) {
        std::cout << "min spatial velocity: " << min_vel << '\n'
              << "max spatial velocity: " << max_vel << '\n'
              << "mean spatial velocity: " << mean_vel << '\n';
    }
    
    // convert to radians:
    double min_lat_deg = lat_deg[0];
    double min_lon_deg = lon_deg[0];
    double max_lat_deg = lat_deg.back();
    double max_long_deg = lon_deg.back();
    double dlat_deg = lat_deg[1]-lat_deg[0];
    double dlon_deg = lon_deg[1]-lon_deg[0];
    double delta_long_deg = max_long_deg - min_lon_deg;
    double delta_lat_deg = max_lat_deg - min_lat_deg;
    
    double mid_lon_deg = 0.5*(lon_deg[0] + lon_deg.back());
    double mid_lat_deg = 0.5*(lat_deg[0] + lat_deg.back());
    
    if (save_mesh) {
        std::vector<vec2> vertices(nlat*nlon);
        for (int i=0; i<nlat*nlon; ++i) {
            double _lon_deg = lon_deg[i%nlon]-mid_lon_deg;
            double _lat_deg = lat_deg[i/nlon];
            vertices[i][0] = Earth_radius*deg2rad(_lon_deg)*cos(deg2rad(_lat_deg));
            vertices[i][1] = Earth_radius*deg2rad(_lat_deg);
        }
        
        if (verbose) {
            double width = Earth_radius*cos(deg2rad(mid_lat_deg))*deg2rad(delta_long_deg);
            double height = Earth_radius*deg2rad(delta_lat_deg);
            std::cout << "region size at center: longitude: " << width/1000. << " km, latitude: " << height/1000. << " km\n";
        }
        
        double vedge = Earth_radius*deg2rad(dlat_deg);
        double hedge_bottom = Earth_radius*cos(deg2rad(min_lat_deg))*deg2rad(dlon_deg);
        double hedge_top = Earth_radius*cos(deg2rad(max_lat_deg))*deg2rad(dlon_deg);
        
        if (verbose) {
            std::cout << "length of vertical edges: " << vedge << '\n';
            std::cout << "length of horizontal edges: from " << hedge_bottom 
                << " to " << hedge_top << ", ratio: " << 100.*hedge_top/hedge_bottom
                << '\n';
        }
        
        vtkStructuredGrid* grid = vtkStructuredGrid::New();
        vtkPoints* points = vtk_utils::make_vtkpoints(vertices);
        grid->SetPoints(points);
        grid->SetDimensions(nlon, nlat, 1);
        
        vtkDataSetWriter* writer = vtkDataSetWriter::New();
        writer->SetInputData(grid);
        writer->SetFileName( (name_out + "-grid.vtk").c_str() );
        writer->SetFileTypeToBinary();
        writer->Write();
        
        grid->Delete();
        writer->Delete();
    }
    
    // convert latitude and longitude from degrees to radians
    lat_rad.resize(lat_deg.size());
    lon_rad.resize(lon_deg.size());
    std::copy(lat_deg.begin(), lat_deg.end(), lat_rad.begin());
    std::copy(lon_deg.begin(), lon_deg.end(), lon_rad.begin());
    std::for_each(lat_rad.begin(), lat_rad.end(), [](double& a){ a=deg2rad(a); });
    std::for_each(lon_rad.begin(), lon_rad.end(), [](double& a){ a=deg2rad(a); });
    
    double err=0;
    for (int i=0; i<100; ++i) {
        err += fabs(deg2rad(lat_deg[i])-lat_rad[i]);
    }
    if (verbose) {
        std::cout << "error=" << err << '\n';
        std::cout << "latitude: ";
        check_delta(lat_deg);
        std::cout << "longitude: ";
        check_delta(lon_deg);
    }
    
    // NEW: oblate Earth model
    oblate ob;
    /* formulae to convert spatial velocity (eastward and northward) in m/s to 
       angular velocity (latitude and longitude) in °/s
        * vy = R(lat_rad)*vlat_rad : vlat_rad: velocity in radians/s.
        * vx = r(lat_rad)*vlon_rad : vlon_rad: velocity in radians/s.
       with R: radius of oblate Earth model in meters, r: radius at given latitude
       Hence: vlat_deg = rad2deg(vy/R(lat_rad)) and vlon_deg = rad2deg(vx/(r(lat_rad)*cos(lat_rad)))
    */
    vel_angular.resize(vel_spatial.size()); 
    spurt::ProgressDisplay progress;
    progress.begin(nlat*nlon, "Computing angular velocity");
    for (int j=0; j<nlat; ++j) {
        for (int i=0; i<nlon; ++i) {
            progress.update(i+j*nlon);
            const vec2& vs = vel_spatial[i+j*nlon];
            vec2& va = vel_angular[i+j*nlon];
            if (vs[0] == _invalid_) {
                va[0] = va[1] = 0.; // zero velocity on land
            }
            else {
                double vlon_rad = vs[0]/ob.parallel_radius(lat_rad[j]);
                double vlat_rad = vs[1]/ob.radius(lat_rad[j]);
                va[0] = rad2deg(vlon_rad);
                va[1] = rad2deg(vlat_rad);
            }
        }
    }
    progress.end();

    // compute vorticity in angular domain (the "spatial" units cancel out)
    std::vector<double> vorticity(vel_spatial.size(), 0); // dvlat/dlon - dvlon/dlat
    double lon_d = lon_rad[1]-lon_rad[0];
    double lat_d = lat_deg[1]-lat_deg[0];
    size_t di = 1;
    size_t dj = nlon;
    progress.begin(nlat*nlon, "Computing vorticity");
    for (int j=0; j<nlat; ++j) {
        for (int i=0; i<nlon; ++i) {
            size_t id = i + j*nlon;
            progress.update(id);
            if (vel_spatial[id][0] == _invalid_) continue;
            double dvlat_dlon = 0;
            double dvlon_dlat = 0;
            bool blocked_left   = i==0      || vel_spatial[id-di][0] == _invalid_;
            bool blocked_right  = i==nlon-1 || vel_spatial[id+di][0] == _invalid_;
            bool blocked_bottom = j==0      || vel_spatial[id-dj][0] == _invalid_;
            bool blocked_top    = j==nlat-1 || vel_spatial[id+dj][0] == _invalid_;
            if ((blocked_left && blocked_right) || (blocked_bottom && blocked_top)) {
                continue;
            }
            if (blocked_left && blocked_bottom) {
                dvlat_dlon = vel_angular[id+di][1] - vel_angular[id][1];
                dvlon_dlat = vel_angular[id+dj][0] - vel_angular[id][0];
            }
            else if (blocked_left && blocked_top) {
                dvlat_dlon = vel_angular[id+di][1] - vel_angular[id][1];
                dvlon_dlat = vel_angular[id][0] - vel_angular[id-dj][0];
            }
            else if (blocked_right && blocked_bottom) {
                dvlat_dlon = vel_angular[id][1] - vel_angular[id-di][1];
                dvlon_dlat = vel_angular[id+dj][0] - vel_angular[id][0];
            }
            else if (blocked_right && blocked_top) {
                dvlat_dlon = vel_angular[id][1] - vel_angular[id-di][1];
                dvlon_dlat = vel_angular[id][0] - vel_angular[id-dj][0];
            }
            else if (blocked_left) {
                dvlat_dlon = vel_angular[id+di][1] - vel_angular[id][1];
                dvlon_dlat = 0.5*(vel_angular[id+dj][0] - vel_angular[id-dj][0]);
            }
            else if (blocked_right) {
                dvlat_dlon = vel_angular[id][1] - vel_angular[id-di][1];
                dvlon_dlat = 0.5*(vel_angular[id+dj][0] - vel_angular[id-dj][0]);
            }
            else if (blocked_bottom) {
                dvlat_dlon = 0.5*(vel_angular[id+di][1] - vel_angular[id-di][1]);
                dvlon_dlat = vel_angular[id+dj][0] - vel_angular[id][0];
            }
            else if (blocked_top) {
                dvlat_dlon = 0.5*(vel_angular[id+di][1] - vel_angular[id-di][1]);
                dvlon_dlat = vel_angular[id][0] - vel_angular[id-dj][0];
            }
            else if (kernel_name == "Scharr" || kernel_name == "scharr") {
                // Scharr operators (see https://en.wikipedia.org/wiki/Sobel_operator#Alternative_operators)
                dvlat_dlon = 
                     -3.*vel_angular[id-di-dj][1] +  3.*vel_angular[id+di+dj][1] 
                    -10.*vel_angular[id-di   ][1] + 10.*vel_angular[id+di   ][1]
                     -3.*vel_angular[id-di+dj][1] +  3.*vel_angular[id+di+dj][1];
                dvlon_dlat = 
                    -3.*vel_angular[id-di-dj][0] - 10.*vel_angular[id-dj   ][0] - 3.*vel_angular[id+di-dj][0]
                    +3.*vel_angular[id-di+dj][0] + 10.*vel_angular[id+dj   ][0] + 3.*vel_angular[id+di+dj][0];
                dvlat_dlon /= 32.;
                dvlon_dlat /= 32.;
            }
            else if (kernel_name == "Sobel" || kernel_name == "sobel") {
                // Sobel operators (see https://en.wikipedia.org/wiki/Sobel_operator)
                dvlat_dlon = 
                     -1.*vel_angular[id-di-dj][1] + 1.*vel_angular[id+di+dj][1] 
                     -2.*vel_angular[id-di   ][1] + 2.*vel_angular[id+di   ][1]
                     -1.*vel_angular[id-di+dj][1] + 1.*vel_angular[id+di+dj][1];
                dvlon_dlat = 
                     -1.*vel_angular[id-di-dj][0] - 2.*vel_angular[id-dj   ][0] - 1.*vel_angular[id+di-dj][0]
                     +1.*vel_angular[id-di+dj][0] + 2.*vel_angular[id+dj   ][0] + 1.*vel_angular[id+di+dj][0];
                dvlat_dlon /= 8.;
                dvlon_dlat /= 8.;
            }
            else if (kernel_name == "finite") {
                dvlat_dlon = vel_angular[id+di][1] - vel_angular[id-di][1];
                dvlon_dlat = vel_angular[id+dj][0] - vel_angular[id-dj][0];
                dvlat_dlon /= 2.;
                dvlon_dlat /= 2.;
            }
            else {
                throw std::runtime_error("Invalid kernel name: " + kernel_name);
            }
            dvlat_dlon /= lon_d;
            dvlon_dlat /= lat_d;
            vorticity[id] = dvlat_dlon-dvlon_dlat;
        }
    }
    progress.end();
    
    std::vector<size_t> sz(3);
    std::vector<double> spc(3);
    std::vector<double> min(3);
    std::vector<int> center(3);
    std::string empty;
    sz[0] = 2;
    sz[1] = nlon;
    sz[2] = nlat;
    spc[0] = AIR_NAN;
    spc[1] = lon_deg[1]-lon_deg[0];
    spc[2] = lat_deg[1]-lat_deg[0];
    min[0] = AIR_NAN;
    min[1] = min_lon_deg;
    min[2] = min_lat_deg;
    center[0] = nrrdCenterUnknown;
    center[1] = nrrdCenterNode;
    center[2] = nrrdCenterNode;
    if (save_uv) {
        spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&vel_spatial[0]), 
            name_out+"-spatial_velocity.nrrd", /*nrrdTypeDouble,*/ sz, spc, min, center, empty);
    }
    spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&vel_angular[0]),
            name_out+"-angular_velocity.nrrd", /*nrrdTypeDouble,*/ sz, spc, min, center, empty);
    spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&vorticity[0]),
            name_out+"-angular_vorticity.nrrd", /*nrrdTypeDouble,*/ 
            std::vector<size_t>(sz.begin()+1, sz.end()), 
            std::vector<double>(spc.begin()+1, spc.end()), 
            std::vector<double>(min.begin()+1, min.end()), 
            std::vector<int>(center.begin()+1, center.end()), empty);
    
    for (size_t i=0; i<names.size(); ++i) {
        spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&variables[i][0]),
            name_out + "-" + names[i], /*nrrdTypeDouble,*/ 
            std::vector<size_t>(sz.begin()+1, sz.end()), 
            std::vector<double>(spc.begin()+1, spc.end()), 
            std::vector<double>(min.begin()+1, min.end()), 
            std::vector<int>(center.begin()+1, center.end()), empty);
    }
    std::cout << "t=" << t << '\n';
            
    return 0;
}
