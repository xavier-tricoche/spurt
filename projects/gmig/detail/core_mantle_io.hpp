template<typename Scalar_, typename Index_>
typename xavier::gmig::core_mantle::Vertex<Scalar_,Index_>::bbox_t
xavier::gmig::core_mantle::
read_text(std::vector<Vertex<Scalar_,Index_> >& data, 
          const std::string& file_name, bool verbose) {
    typedef Scalar_                    scalar_t;
    typedef Index_                     index_t;
    typedef Vertex<scalar_t,index_t>   vertex_t;
    typedef typename vertex_t::pos_t   pos_t;
    typedef typename vertex_t::bbox_t  bbox_t;
    
    data.clear();
    bbox_t _bounds;
    std::fstream file(file_name.c_str(), std::ios::in);

    if (file.fail()) {
        throw std::runtime_error("unable to open " + file_name);
    }
    
    while (!file.eof()) {
        if (file.peek() == '#') { // skip comment line
            std::string line;
            std::getline(file, line);
            continue;
        }
        
        vertex_t _vert;
        file >> _vert.layer_id
             >> _vert.position[0] 
             >> _vert.position[1]
             >> _vert.position[2] 
             >> _vert.type
             >> _vert.lat_sector
             >> _vert.long_sector;
        if (!_vert.valid()) {
            // std::cerr << "invalid vertex information on line " << n << ": " << _vert << '\n';
            if (file.fail() || file.bad()) {
                // std::cerr << "reading error. stopping file import on line " << n << ".\n";
                break;
            }
            else continue;
        }
        
        data.push_back(_vert);
        _bounds.add(_vert.position);
    }
    file.close();
    
    return _bounds;
}
