// stl
#include <sstream>
// xavier
#include <image/nrrd_wrapper.hpp>
#include "dataIO.hpp"

nvis::bbox2 xavier::gmig::
read_nrrd(std::vector<nvis::vec2>& points, 
          std::vector<nvis::vec1>& times,
          std::vector<nvis::vec1>& weights,
          double& velocity,
          nvis::vec3& source,
          std::string& kernel_name,
          const std::string& filename,
          bool verbose = false) {
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        std::cerr << "read_data_from_nrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    std::vector<double> data;
    xavier::to_vector(data, nin);
    size_t N = nin->axis[0].size;
    bool has_weights = (N == 4);
    size_t nb_pts = data.size()/N;
    points.resize(nb_pts);
    times.resize(nb_pts);
    if (has_weights) weights.resize(nb_pts);
    if (verbose && has_weights) 
        std::cout << "input data contains precomputed weights\n";
    nvis::bbox2 _bounds;
    for (size_t i=0 ; i<nb_pts ; ++i) {
        points[i][0] = data[N*i  ];
        points[i][1] = data[N*i+1];
        times[i][0]  = data[N*i+2];
        if (has_weights) weights[i][0] = data[N*i+3];
        _bounds.add(points[i]);
    }
    nrrdNuke(nin);
    return _bounds;
}

void xavier::gmig::
save_rbf(const std::vector<nvis::vec2>& points,
         const std::vector<nvis::vec1>& times,
         const std::vector<nvis::vec1>& weights,
         double velocity,
         const nvis::vec2& source,
         const std::string& kernel_name,
         const std::string& filename,
         bool verbose = false) {
    size_t npts = points.size();
    double *data = (double*)calloc(npts*4, sizeof(double));
    for (size_t i=0 ; i<npts ; ++i) {
        data[4*i  ] = points[i][0];
        data[4*i+1] = points[i][1];
        data[4*i+2] = values[i][0];
        data[4*i+3] = weights[i][0];
    }
    xavier::nrrd_params<double, 2> params;
    params.sizes()[0] = 4;
    params.sizes()[1] = npts;
    std::ostringstream os;
    os << "RBF reconstruction data for kernel " << kernel_name
       << " with source at " << source << " and mean isotropic velocity "
       << velocity;
    params.description() = os.str();
    params.labels()[0] = "x_rec;y_rec;time;weight";
    params.labels()[1] = "RBF centers";
    xavier::writeNrrd(data, filename, params);
    if (verbose) std::cout << filename << " has been exported\n";
}