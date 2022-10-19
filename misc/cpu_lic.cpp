#include <iostream>
#include <string>
#include <image/nrrd_wrapper.hpp>
#include <math/fixed_vector.hpp>
#include <vis/integral_curve.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef integral_curve<2> sl_type;

template< typename T >
class image {
public:
    typedef T         value_type;
    typedef T*        iterator_type;
    typedef const T*  const_iterator_type;
    
    image(size_t w, size_t h) : _w(w), _h(h), _size(w*h) {
        _ptr = (value_type*)calloc(_w*_h, sizeof(value_type));
    }
    
    bool valid(int i, int j) const {
        return (i<0 || i>=_w || j<0 || j>=_h);
    }
    
    value_type& operator()(int i, int j) {
        if (!valid(i,j)) 
            throw std::runtime_error("index out of bounds");
        return _ptr[i+j*_w];
    }
    const value_type& operator()(int i, int j) const {
        if (!valid(i,j))
            throw std::runtime_error("index out of bounds");
        return _ptr[i+j*_w];
    }
    
    iterator_type begin() {
        return iterator_type(_ptr);
    }
    
    iterator_type end() {
        return iterator_type(_ptr + _size)
    }

    const_iterator_type begin() const {
        return const_iterator_type(_ptr);
    }
    
    const_iterator_type end() const {
        return const_iterator_type(_ptr + _size)
    }
    
    size_t size() const { return _size; }
    
    ivec2 dimensions() const {
        return ivec2(_w, _h);
    }
    
private:
    size_t      _w, _h, _size;
    value_type* _ptr;
};


template< typename T1, typename T2 >
class grid {
public:
    typedef T1                                   coord_type;
    typedef typename fixed_vector<T1, 2>   pos_type;
    typedef ivec2                          index_type;
    typedef bounding_box<pos_type>         bounds_type;
    typedef T2                                   value_type;
    typedef array<T2>                            array_type;

public:
    grid(size_t w, size_t h, const box_type& bounds)
        : _array(w, h), _res(w, h), _bounds(bounds) {
        _spacing = _bounds.size() / pos_type(_res - index_type(1,1));
    }
    
    size_t size() const {
        return _array.size();
    }
    
    const bounds_type& bounds() const {
        return _bounds;
    }
    
    const array_type& raster() const {
        return _array;
    }
    
    const pos_type& spacing() const {
        return _spacing;
    }
    
    const index_type& dimensions() const {
        return _res;
    }
    
    const value_type& operator()(const index_type& id) const {
        return _array(id);
    }

    value_type& operator()(const index_type& id) {
        return _array(id);
    }
    
    const index_type& index(const pos_type& p) const {
        return index_type(p-_bounds.min()/_spacing);
    }
    
    value_type interpolate(const pos_type& p) const {
        if (!inside(p)) throw std::runtime_error("position out of bounds");
        pos_type local = (p - _bounds.min())/_spacing();
        index_type id(local);
        local -= pos_type(id); 
        double u = local[0];
        double v = local[1];
        return (1.-u)*(1.-v)*_array(id) + 
               u*(1.-v)*_array(id + index_type(1,0)) +
               u*v*_array(id + index_type(1,1)) +
               (1.-u)*v*_array(id + index_type(0,1));
    }
    
    bool inside(const pos_type& p) const {
        return _bounds.inside(p);
    }
    
private:
    size_t       _size;
    index_type   _res;
    array_type   _data;
    bounds_type  _bounds;
    pos_type     _spacing;
};

void make_noise(image<float>& noise) {
    image<float>::index_type dims = noise.dimensions();
    srand48(time(0));
    
    #pragma omp parallel for
    for (size_t n=0 ; n<noise.size() ; ++n) {
        int i = n % dims[0];
        int j = n / dims[0];
        noise(i,j) = drand48();
    }
}

std::string me;
void usage(const std::string& msg="") {
    if (!msg.empty()) 
        std::cerr << "\nERROR: " << msg << '\n';
    std::cout << '\n'
              << "DESCRIPTION: compute LIC image on the CPU.\n"
              << '\n'
              << "USAGE: " << me << " [parameters] [options]\n"
              << '\n'
              << "PARAMETERS:\n"
              << " -i | --input <string>      Input filename\n"
              << " -o | --output <string>     Output filename\n"
              << '\n'                         
              << "OPTIONS:\n"
              << " -h | --help                Print this information\n"
              << " -l | --length <int>        Integration length\n"
              << " -r | --resolution <int>x2  Resolution\n"
              << " -e | --epsilon <float>     Integration precision\n"
              << " -v | --verbose             Turn on verbose mode\n"
              << '\n';
              
     exit(!msg.empty());
}

int main(int argc, char* argv[]) {
    
    typedef float                       value_type;
    typedef fvec2                 vec_type;
    typedef grid<vec_type, vec_type>    grid_type;
    typedef grid_type::bounds_type      bounds_type;
    typedef grid_type::array_type       array_type;
    
    int w=0, h=0, len=100;
    double eps=1.0e-6;
    std::string in_name="", out_name="";
    bool verbose = false;
    
    me = argv[0];
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-v" || arg == "--verbose") verbose = true;
        else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) usage("missing input filename");
            in_name = argv[++i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) usage("missing output filename");
            out_name = argv[++i];
        }
        else if (arg == "-l" || arg == "--length") {
            if (i == argc-1) usage("missing integration length");
            len = atoi(argv[++i]);
        }
        else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) usage("missing resolution information");
            w = atoi(argv[++i]);
            h = atoi(argv[++i]);
        }
        else if (arg == "-e" || arg == "--epsilon") {
            if (i == argc-1) usage("missing integration precision");
            eps = atof(argv[++i]);
        }
    }
    
    if (in_name.empty()) usage("no input file provided");
    if (out_name.empty()) usage("no output file provided");
    
    Nrrd* nin = spurt::readNrrd(in_name.c_str());
    bounds_type bounds = spurt::compute_bounds(nin);
    std::vector<value_type> data;
    spurt::to_vector(data, nin);
    if (w*h == 0) {
        w = nin->axis[1].size;
        h = nin->axis[2].size;
    }
    
    grid_type rhs(nin->axis[1].size, nin->axis[2].size, bounds);
    #pragma openmp parallel for
    for (int n=0 ; n<rhs.size() ; ++n) {
        int i = n % rhs.dimensions()[0];
        int j = n / rhs.dimensions()[1];
        rhs(i,j) = vec_type(data[2*n], data[2*n+1]);
    }
    
    
    
}

