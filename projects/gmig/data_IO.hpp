#ifndef __XAVIER_COLLAB_MAARTEN_DATA_IO_HPP__
#define __XAVIER_COLLAB_MAARTEN_DATA_IO_HPP__

// STL
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt { namespace maarten {

// custom reader for Hui Huang's files
// format:
// source_lon; source_lat; receiver_lon; receiver_lat; travel_time; \
//             <optional: source-receiver distance>; receiver index.
template<typename Time_=double, typename Dist_=double>
nvis::bbox2 read_text(std::vector<nvis::vec2>& points, 
                      std::vector<Time_>& times,
                      std::vector<Dist_>& distances,
                      const std::string& file_name,
                      bool verbose=false) {
    static const double invalid = std::numeric_limits<double>::max();
    points.clear();
    times.clear();
    distances.clear();
    nvis::bbox2 _bounds;
    double min=invalid, max=-invalid;
    std::fstream file(file_name.c_str(), std::ios::in);

    if (file.fail()) {
        throw std::runtime_error("unable to open " + file_name);
    }

    // check if file contains information about source-receiver distance
    bool has_7_terms = false;
    {
        std::string first_line;
        std::getline(file, first_line);
        std::istringstream iss(first_line);
        size_t nwords;
        for (nwords=0 ; true ; ++nwords) {
            std::string word;
            iss >> word;
            if (word.empty()) break;
        }
        has_7_terms = ( nwords == 7 );
        if (verbose) {
            std::cout << "data file has " << nwords << " terms per row\n";
        }
        file.seekg(0);
    }

    while (!file.eof()) {
        double sx, sy, x, y, t, d, sentinel=invalid;
        file >> sx >> sy >> x >> y >> t;
        if (points.empty()) {
            // add data corresponding to the source
            points.push_back(nvis::vec2(sx, sy));
            times.push_back(Time_(0));
            if (has_7_terms) distances.push_back(Dist_(0));
        }
        if (has_7_terms) file >> d;
        file >> sentinel;
        if (sentinel == invalid) break;
        points.push_back(nvis::vec2(x,y));
        times.push_back(Time_(t));
        if (has_7_terms) distances.push_back(Dist_(d));
        min = std::min(t, min);
        max = std::max(t, max);
        _bounds.add(points.back());
    }
    file.close();

    return _bounds;
}
    
} // maarten
} // spurt

#endif