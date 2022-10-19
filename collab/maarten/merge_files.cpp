#include <fstream>
#include <iostream>
#include <string>
#include <math/fixed_vector.hpp>
#include <misc/option_parse.hpp>
#include <image/nrrd_wrapper.hpp>

// custom reader for travel time data files
// format:
// source_lon; source_lat; receiver_lon; receiver_lat; travel_time; receiver index
void read_file(std::vector<nvis::fvec4>& points, std::vector<float>& times,
               const std::string& file_name) {
    static const float invalid = std::numeric_limits<float>::max();
    std::fstream file(file_name.c_str(), std::ios::in);
    while (!file.eof()) {
        float sx, sy, rx, ry, t, sentinel=invalid;
        file >> sx >> sy >> rx >> ry >> t >> sentinel;
        if (sentinel == invalid) break;
        points.push_back(nvis::fvec4(sx, sy, rx, ry));
        times.push_back(t);
    }
    file.close();
}


int main(int argc, char* argv[]) {
    namespace bpo = boost::program_options;
    namespace clt = xavier::cmdline_tools;

    bpo::options_description desc(
        "Merge a set of travel time data files into a single NRRD file");

    std::string input_name, output_name, path="";
    bool verbose;
    desc.add_options()
    ("help,h",                                                 "Display this message")
    ("input,i",   clt::value(&input_name, false)->required(),  "File containing list of input file names")
    ("output,o",  clt::value(&output_name, false)->required(), "Output file name")
    ("path,p",    clt::value(&path),                           "Path to be prepended to all file names")
    ("verbose,v", bpo::bool_switch(&verbose),                  "Turn ON verbose mode");

    // desc.add_options()
    //     ("help,h",                                                 "Display this message")
    //     ("input,i",   bpo::value<std::string>(&input_name)->required(),  "File containing list of input file names")
    //     ("output,o",  bpo::value<std::string>(&output_name)->required(), "Output file name")
    //     ("path,p",    bpo::value<std::string>(&path),                           "Path to be prepended to all file names")
    //     ("verbose,v", bpo::bool_switch(&verbose),                  "Turn ON verbose mode");

    // bpo::variables_map vm;
    // try {
    //     bpo::store(bpo::parse_command_line(argc, argv, desc), vm);
    //     if (vm.count("help")) {
    //         clt::usage(desc, argv[0]);
    //     }
    //     bpo::notify(vm);
    // } catch(bpo::error& e) {
    //     clt::error(desc, e.what());
    // } catch(std::runtime_error& e) {
    //     clt::error(desc, e.what());
    // }

    clt::parse_command_line(argc, argv, desc);
    if (!path.empty()) {
        if (*path.rbegin() != '/') path.push_back('/');
        output_name = path + output_name;
    }
    std::vector<std::string> input_files;
    std::fstream _file(input_name.c_str(), std::ios::in);
    if (_file.is_open()) {
        std::string name = "invalid!";
        while (!_file.eof()) {
            _file >> name;
            if (name == "invalid!") break;
            input_files.push_back(path + name);
        }
        _file.close();
    }
    else {
        std::cerr << "ERROR: unable to open " + input_name << '\n';
        std::cerr << desc << '\n';
    }
    if (verbose) {
        std::cout << "there were " << input_files.size() << " files in input:\n";
        std::copy(input_files.begin(), input_files.end(),
                  std::ostream_iterator<std::string>(std::cout, "\n"));
    }
    std::vector<nvis::fvec4> points;
    std::vector<float> values;
    for (size_t i=0 ; i<input_files.size() ; ++i) {
        read_file(points, values, input_files[i]);
    }
    if (verbose) {
        std::cout << "there are a total of " << points.size() << " data points\n";
    }

    float* data = (float*)calloc(5*points.size(), sizeof(float));
    for (size_t i=0 ; i<points.size() ; ++i) {
        data[5*i  ] = points[i][0];
        data[5*i+1] = points[i][1];
        data[5*i+2] = points[i][2];
        data[5*i+3] = points[i][3];
        data[5*i+4] = values[i];
    }
    xavier::nrrd_params<float, 2> params;
    params.sizes()[0] = 5;
    params.sizes()[1] = points.size();
    xavier::writeNrrd(data, output_name, params);
    std::cout << output_name << " has been exported\n";
    return 0;
}