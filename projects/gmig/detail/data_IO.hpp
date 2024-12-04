template<typename Scalar_>
nvis::bbox2 xavier::gmig::traveltime::
read_text(travel_time_data<Scalar_>& data,
          const std::string& file_name,
          bool verbose) {
    typedef Scalar_ scalar_type;
    
    data.clear();
    nvis::bbox2 _bounds;
    double min=invalid_double, max=-invalid_double;
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
        double sx, sy, x, y, t, d, sentinel=invalid_double;
        file >> sx >> sy >> x >> y >> t;
        if (data.receivers.empty()) {
            // add data corresponding to the source
            data.receivers.push_back(nvis::vec2(sx, sy));
            data.times.push_back(scalar_type(0));
            if (has_7_terms) data.distances.push_back(scalar_type(0));
            data.source = nvis::vec2(sx, sy);
        }
        if (has_7_terms) file >> d;
        file >> sentinel;
        if (sentinel == invalid_double) break;
        data.receivers.push_back(nvis::vec2(x,y));
        data.times.push_back(scalar_type(t));
        if (has_7_terms) data.distances.push_back(scalar_type(d));
        min = std::min(t, min);
        max = std::max(t, max);
        _bounds.add(data.receivers.back());
    }
    file.close();
    if (!data.sanity_check()) {
        throw std::runtime_error("read_txt: \
        Invalid data travel time data object");
    }
    data.compute_velocity();

    return _bounds;
}

template<typename Scalar_>
void xavier::gmig::traveltime::
save_rbf(const travel_time_data<Scalar_>& data,
         const std::string& filename,
         bool verbose) {
    
    if (!data.sanity_check()) {
        throw std::runtime_error("save_rbf: invalid travel time data");
    }
    size_t data_size = 3;
    bool has_distances = !data.distances.empty();
    bool has_weights = !data.weights.empty();
    bool has_kernel = (data.kernel_name != "none");
    if (!has_kernel && has_weights) 
        throw std::runtime_error(
            "save_rbf: Missing RBF kernel name for available weights"
        );
    if (has_kernel && !has_weights)
        std::cerr << "Warning: save_rbf: "
        << "Missing weights for available RBF kernel name\n";
    if (has_distances) ++data_size;
    if (has_weights) ++data_size;
    
    size_t npts = data.receivers.size();
    double *raw = (double*)calloc(npts*data_size, sizeof(double));
    size_t c = 0;
    for (size_t i=0 ; i<npts ; ++i) {
        raw[c++] = data.receivers[i][0];
        raw[c++] = data.receivers[i][1];
        raw[c++] = data.times[i];
        if (has_distances) raw[c++] = data.distances[i];
        if (has_weights) raw[c++] = data.weights[i];
    }
    xavier::nrrd_params<double, 2> params;
    params.sizes()[0] = data_size;
    params.sizes()[1] = npts;
    params.comments().push_back("Travel time dataset");
    std::ostringstream os;
    os << "source=" << data.source;
    params.comments().push_back(os.str());
    if (has_kernel) {
        os.clear(); os.str("");
        os << "kernel=" << data.kernel_name;
        params.comments().push_back(os.str());
        os.clear(); os.str("");
        if (data.radius > 0) {
            os << "radius=" << data.radius;
            params.comments().push_back(os.str());
            os.clear(); os.str("");
        }
        os << "velocity=" << data.velocity;
        params.comments().push_back(os.str());
    }
    os.clear(); os.str("");
    os << "x_rec;y_rec;time";
    if (has_distances) os << ";distances";
    if (has_weights) os << ";weights";
    params.labels()[0] = os.str();
    params.labels()[1] = "data_points";
    xavier::writeNrrd(raw, filename, params);
    if (verbose) std::cout << filename << " has been exported\n";
}

template<typename Scalar_>
nvis::bbox2 xavier::gmig::traveltime::
read_nrrd(travel_time_data<Scalar_>& data,
          const std::string& filename,
          bool verbose) {
    typedef Scalar_ scalar_type;
              
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL) || nin->dim != 2) {
        throw std::runtime_error(error_msg("read_nrrd"));
    }
    // check axis labels
    bool has_weights = false;
    bool has_distances = false;
    char **tmp_char;
    nrrdAxisInfoGet_nva(nin, nrrdAxisInfoLabel, tmp_char);
    std::string label1(tmp_char[0]);
    std::string label2(tmp_char[1]);
    // 4 cases for 1st label:
    // 1. "x_rec;y_rec;time": 16 characters
    // 2. "x_rec;y_rec;time;distances": 26 characters 
    // 3. "x_rec;y_rec;time;weights": 24 characters 
    // 4. "x_rec;y_rec;time;distances;weights": 34 characters
    if (label1 != "x_rec;y_rec;time" &&
        label1 != "x_rec;y_rec;time;distances" &&
        label1 != "x_rec;y_rec;time;weights" &&
        label1 != "x_rec;y_rec;time;distances;weights") {
        throw std::runtime_error
            ("read_nrrd: invalid 1st label: " + label1);
    }
    if (label1.size() == 26) has_distances = true;
    else if (label1.size() == 24) has_weights = true;
    else if (label1.size() == 34) has_distances = has_weights = true;
    
    // check comments
    if (!nin->cmtArr->len || 
        std::string(nin->cmt[0]) != "Travel time dataset") {
        throw std::runtime_error("read_nrrd: invalid comments");
    }
    std::vector<double> raw;
    xavier::to_vector(raw, nin);
    size_t N = nin->axis[0].size;
    size_t alt_N = 3;
    if (has_distances) ++alt_N;
    if (has_weights) ++alt_N;
    if (N != alt_N) {
        throw std::runtime_error("read_nrrd: Inconsistent dimensions");
    }
    size_t nb_pts = data.size()/N;
    data.receivers.resize(nb_pts);
    data.times.resize(nb_pts);
    if (has_distances) data.distances.resize(nb_pts);
    if (has_weights) data.weights.resize(nb_pts);
    if (verbose && has_weights) 
        std::cout << "input data contains precomputed weights\n";
    nvis::bbox2 _bounds;
    size_t c=0;
    for (size_t i=0 ; i<nb_pts ; ++i) {
        data.receivers[i][0] = raw[c++];
        data.receivers[i][1] = raw[c++];
        _bounds.add(data.receivers[i]);
        data.times[i] = scalar_type(raw[c++]);
        if (has_distances) data.distances[i] = scalar_type(raw[c++]);
        if (has_weights) data.weights[i]= scalar_type(raw[c++]);
    }
    bool has_kernel_name = false;
    bool has_source_location = false;
    std::istringstream iss (std::string(nin->cmt[1]));
    std::string tmp;
    // read the remainder of the comments
    for (size_t i=1; i<nin->cmtArr->len ; ++i) {
        iss.clear();
        iss.str(std::string(nin->cmt[i]));
        iss >> tmp;
        if (iss.fail()) throw std::runtime_error
            ("read_nrrd: invalid comment item");
        if (tmp == "source=") {
            iss >> data.source;
            if (iss.fail()) throw std::runtime_error
                ("read_nrrd: invalid source information" + iss.str());
            else has_source_location = true;
        }
        if (tmp == "velocity=") {
            iss >> data.velocity;
            if (iss.fail()) throw std::runtime_error
                ("read_nrrd: invalid velocity data: " + iss.str());
        }
        else if (tmp == "radius=") {
            iss >> data.radius;
            if (iss.fail()) throw std::runtime_error
                ("read_nrrd: invalid radius data: " + iss.str());
        }
        else if (tmp.substr(0, 7) == "kernel=") {
            data.kernel_name = tmp.substr(7);
        }
        else if (verbose){
            std::cout << "read_nrrd: Warning: ignoring unrecognized comment"
                << std::endl;
        }
    }
    nrrdNuke(nin);
    
    if (!data.sanity_check()) {
        throw std::runtime_error("read_nrrd: \
        Invalid data travel time data object");
    }
    
    data.compute_velocity();
    
    return _bounds;
}
