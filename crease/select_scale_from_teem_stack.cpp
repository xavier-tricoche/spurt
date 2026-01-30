// std
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <tuple>
#include <regex>
#include <mutex>
#include <thread>
// Boost filesystem (c++17 compiler flaky with std::filesystem )
#include <boost/filesystem.hpp>
// Teem
#include <teem/nrrd.h>
// spurt
#include <math/types.hpp>
#include <math/stat.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <data/image.hpp>
#include <misc/cxxopts.hpp>
// TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>

using namespace spurt;
typedef double scalar_type;
typedef long size_type;
typedef small_vector<scalar_type, 3> pos_type;
typedef small_vector<size_type, 3> coord_type;
typedef small_vector<scalar_type, 3> vec_type;
typedef small_matrix<scalar_type, 3, 3> mat_type;

typedef raster_data<size_type, scalar_type, 3, scalar_type> scalar_raster;
typedef raster_data<size_type, scalar_type, 3, vec_type> vector_raster;
typedef raster_data<size_type, scalar_type, 3, mat_type> matrix_raster;

typedef image<size_type, scalar_type, 3, scalar_type, kernels::MitchellNetravaliBC> smooth_image;

namespace fs = boost::filesystem;

Nrrd* import_nrrd(const std::string& filename)
{
    Nrrd* nin = nrrdNew();
    std::cout << "importing " << filename << '\n';
    if (nrrdLoad(nin, filename.c_str(), NULL))
    {
        char* err = biffGetDone(NRRD);
        std::cerr << "Thomas Schultz's ridge method: " << err << std::endl;
        exit(-1);
    }
    return nin;
}

scalar_raster import_nrrd_as_raster(const std::string& filename) 
{
    Nrrd *nin = import_nrrd(filename);
    return nrrd_utils::to_raster<size_type, scalar_type, 3, scalar_type>(nin, true);
}

scalar_type default_lmin_treshold(const scalar_raster& values) {
    scalar_type mode = 
        spurt::mode(values.begin(), values.end(), 1024, 
                    std::numeric_limits<scalar_type>::min(), scalar_type(0));
    std::cout << "mode is " << mode << '\n';
    return mode;
}

struct scale_info {
    double scale;
    std::string strength, gradient, hessian;
};

std::vector<scale_info> import_scales(const std::string& basename, 
                                      const std::string& str_name, 
                                      const std::string& gv_name, 
                                      const std::string& hess_name) {
    std::vector<scale_info> scales;

    std::cout << "basename is " << basename << '\n';

    fs::path apath(basename);
    fs::path directory = apath.parent_path();

    std::cout << "directory is " << directory << '\n';
    fs::path search_path = apath.filename().string();
    std::string search_str = search_path.string();
    std::regex base_pattern(search_str);
    std::cout << "search str is " << search_str << '\n';

    std::string scale_pattern = "([0-9]+[.][0-9]+)";
    std::cout << "scale pattern is " << scale_pattern << '\n';

    std::regex str_pattern(str_name + "_scale_" + scale_pattern + "\\.nrrd");
    std::cout << "str pattern is " << str_name + "_scale_" + scale_pattern + "\\.nrrd" << '\n';
    std::regex gv_pattern(gv_name + "_scale_" + scale_pattern + "\\.nrrd");
    std::cout << "gv pattern is " << gv_name + "_scale_" + scale_pattern + "\\.nrrd" << '\n';
    std::regex hess_pattern(hess_name + "_scale_" + scale_pattern + "\\.nrrd");
    std::cout << "hess pattern is " << hess_name + "_scale_" + scale_pattern + "\\.nrrd" << '\n';

    std::map<float, std::string> scales_to_str;
    std::map<float, std::string> scales_to_gv;
    std::map<float, std::string> scales_to_hess;

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (fs::is_regular_file(entry)) {
            std::string filename = entry.path().filename().string();
            std::cout << "filename is " << filename << '\n';
            if (std::regex_search(filename, base_pattern)) {
                std::smatch match;
                if (std::regex_search(filename, match, str_pattern)) {
                    std::cout << "found scale " << match[1].str() << " in strength " << filename << '\n';
                    scales_to_str[std::stof(match[1].str())] = entry.path().string();
                }
                else if (std::regex_search(filename, match, gv_pattern)) {
                    std::cout << "found scale " << match[1].str() << " in gradient " << filename << '\n';
                    scales_to_gv[std::stof(match[1].str())] = entry.path().string();
                }
                else if (std::regex_search(filename, match, hess_pattern)) {
                    std::cout << "found scale " << match[1].str() << " in hessian " << filename << '\n';
                    scales_to_hess[std::stof(match[1].str())] = entry.path().string();
                }
            }
        }
    }
    for (auto it = scales_to_str.begin(); it != scales_to_str.end(); ++it) {
        auto s = it->first; 
        auto str = it->second;
        auto gv = scales_to_gv[s];
        auto hess = scales_to_hess[s];

        if (gv == "" || hess == "") {
            std::cout << "skipping scale " << s << " because gradient or hessian is missing\n";
            continue;
        }
        scales.push_back(scale_info());
        scales.back().scale = s;
        scales.back().strength = str;
        scales.back().gradient = gv;
        scales.back().hessian = hess;
    }
    return scales;
}

void update_rasters(scalar_raster& best_strength, vector_raster& best_gradient, 
                    matrix_raster& best_hessian, scalar_raster& best_scale, 
                    const scalar_raster& strength, const vector_raster& gradient, 
                    const matrix_raster& hessian, scalar_type scale) {
    ProgressDisplay progress;
    std::atomic<size_t> tbb_progress_counter;
    std::mutex update_progress_mutex;

    std::cout << std::endl;
    progress.begin(best_strength.size(), "Updating rasters at scale " + std::to_string(scale));
    tbb::parallel_for(tbb::blocked_range<size_t>(0,best_strength.size()),
                       [&](tbb::blocked_range<size_t> r) {
    for (size_t n=r.begin(); n!=r.end(); ++n)
    {
        tbb_progress_counter++;

        if (update_progress_mutex.try_lock() ){
            progress.update(tbb_progress_counter);
            update_progress_mutex.unlock();
        }
        if (strength[n] < best_strength[n]) {
            best_strength[n] = strength[n];
            best_gradient[n] = gradient[n];
            best_hessian[n] = hessian[n];
            best_scale[n] = scale;  
        }
    }});
    progress.end();
}

int main(int argc, const char* argv[]) 
{
    std::string input_name, output_name, path;

    cxxopts::Options options("filter_crease", "Manipulate crease surface");
    options.add_options()
        ("i,input", "Input filenames regex with %s for what %f for scale", cxxopts::value<std::string>())
        ("o,output", "Output basename", cxxopts::value<std::string>())
        ("p,path", "Path to prepend to input and outputfiles", cxxopts::value<std::string>())
        ("gradname", "Gradient value name", cxxopts::value<std::string>()->default_value("gv"))
        ("strname", "Ridge strength name", cxxopts::value<std::string>()->default_value("heval2"))
        ("hessname", "Hessian name", cxxopts::value<std::string>()->default_value("hess"))
        ("v,verbose", "Verbose output", cxxopts::value<bool>())
        ("h,help", "Print usage information");
    
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("input") || !result.count("output")) {
        std::cout << options.help() << '\n';
        exit(0);
    }

    bool verbose = result.count("verbose");
    std::string input = result["input"].as<std::string>();
    std::string output = result["output"].as<std::string>();
    std::string str_name = result["strname"].as<std::string>();
    std::string gv_name = result["gradname"].as<std::string>();
    std::string hess_name = result["hessname"].as<std::string>();
    if (result.count("path")) { 
        std::string path = result["path"].as<std::string>();
        fs::path input_path(path);
        fs::path output_path(path);
        input_path /= input;
        output_path /= output; 
        input = input_path.string();
        output = output_path.string(); 
        std::cout << "input=" << input << " output=" << output << '\n';
    }

    std::vector<scale_info> scales = import_scales(input, str_name, gv_name, hess_name); 

    if (scales.empty()) {
        std::cerr << "No scales found\n";
        exit(-1);
    }

    // import a reference image for dimensions
    scalar_raster ref = import_nrrd_as_raster(scales[0].strength);
    auto grid = ref.grid();
    scalar_raster best_strength = scalar_raster(grid, 1.); // 1 is an invalid ridge strength value
    vector_raster best_gradient = vector_raster(grid);
    matrix_raster best_hessian = matrix_raster(grid);
    scalar_raster best_scale = scalar_raster(grid, 0.);
    
    std::cout << "found " << scales.size() << " scales\n";
    for (auto& info : scales) {
        scalar_type scale = info.scale;
        scalar_raster strength = nrrd_utils::to_raster<size_type, scalar_type, 3, scalar_type>(import_nrrd(info.strength), true);
        vector_raster gradient = nrrd_utils::to_raster<size_type, scalar_type, 3, vec_type>(import_nrrd(info.gradient), false);
        matrix_raster hessian = nrrd_utils::to_raster<size_type, scalar_type, 3, mat_type>(import_nrrd(info.hessian), false);
        update_rasters(best_strength, best_gradient, best_hessian, best_scale, strength, gradient, hessian, scale);
    }

    save_as_nrrd(output + "scalespace_strength.nrrd", best_strength);
    save_as_nrrd(output + "scalespace_gradient.nrrd", best_gradient);
    save_as_nrrd(output + "scalespace_hessian.nrrd", best_hessian);
    save_as_nrrd(output + "scalespace_scale.nrrd", best_scale);

    return 0;
}