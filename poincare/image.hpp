#ifndef __IMAGE_HPP__
#define __IMAGE_HPP__

#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <iostream>

struct SubImage {

    SubImage(const nvis::bbox2& inner,
             const nvis::bbox2& outer, int M, int N)
        : _bbox(inner) {
        
        double hx = (outer.max()[0] - outer.min()[0]) / (double)M;
        double hy = (outer.max()[1] - outer.min()[1]) / (double)N;
        
        _ibbox.min()[0] = std::max(0, std::min(M - 1, (int)floor(_bbox.min()[0] / hx)));
        _ibbox.min()[1] = std::max(0, std::min(N - 1, (int)floor(_bbox.min()[1] / hy)));
        _ibbox.max()[0] = std::max(0, std::min(M - 1, (int)floor(_bbox.max()[0] / hx)));
        _ibbox.max()[1] = std::max(0, std::min(N - 1, (int)floor(_bbox.max()[1] / hy)));
        
        std::cout << "integral bounding box is " << _ibbox << std::endl;
    }
    
    bool inside(unsigned int u, unsigned int v) const {
        return _ibbox.inside(nvis::ivec2(u, v));
    }
    
    bool inside(const nvis::vec2& p) const {
        return _bbox.inside(p);
    }
    
    nvis::bounding_box<nvis::vec2> _bbox;
    nvis::bounding_box<nvis::ivec2> _ibbox;
    
};

template< typename T >
class Image {
public:
    void map_coord(int& i, int& j) const {
        if (periodic[0]) {
            if (i < 0) {
                i += m;
            } else if (i > m - 1) {
                i -= m;
            }
        }
        if (periodic[1]) {
            if (j < 0) {
                j += n;
            } else if (j > n - 1) {
                j -= n;
            }
        }
        if (!(i >= 0 && i < m && j >= 0 && j < n)) {
            std::cout << "DAMN!" << std::endl;
        }
    }
    
    Image(const Image& image)
        : m(image.m), n(image.n), data(image.data), shift(0) {}
        
    Image(Nrrd* nrrd, bool is_periodic[2])
        : m(nrrd->axis[0].size), n(nrrd->axis[1].size),
          data((T*)nrrd->data), shift(0) {
        periodic[0] = is_periodic[0];
        periodic[1] = is_periodic[1];
    }
    
    Image(T* _data, unsigned int _m, unsigned int _n, bool is_periodic[2])
        : m(_m), n(_n), data(_data), shift(0) {
        periodic[0] = is_periodic[0];
        periodic[1] = is_periodic[1];
    }
    
    const T& operator()(int i, int j) const {
        map_coord(i, j);
        return data[i+j*m+shift];
    }
    T& operator()(int i, int j) {
        map_coord(i, j);
        return data[i+j*m+shift];
    }
    
    void set_slice(unsigned int s) {
        shift = s * m * n;
    }
    
    const T* get_data() {
        return &data[shift*m*n];
    }
    
    T* data;
    unsigned int m, n;
    bool periodic[2];
    unsigned int shift;
};


template< typename T >
class ScalarImage : public Image<T> {
public:
    ScalarImage(Nrrd* nrrd, bool is_periodic[2])
        : Image<T>(nrrd, is_periodic) {
        std::cerr << " m = " << this->m << " x n = " << this->n << std::endl;
    }
    
    ScalarImage(T* _data, unsigned int _m, unsigned int _n, bool is_periodic[2])
        : Image<T>(_data, _m, _n, is_periodic) {}
        
    nvis::vec2 gradient(int i, int j) const {
        const ScalarImage<T>& self = *this;
        return 0.5*nvis::vec2(self(i + 1, j) - self(i - 1, j),
                              self(i, j + 1) - self(i, j - 1));
    }
    nvis::vec3 hessian(int i, int j) const {
        nvis::vec2 dgdx = 0.5 * (gradient(i + 1, j) - gradient(i - 1, j));
        nvis::vec2 dgdy = 0.5 * (gradient(i, j + 1) - gradient(i, j - 1));
        return nvis::vec3(dgdx[0], dgdx[1], dgdy[1]);
    }
};

template< typename T >
void xy2ij(int& i, int& j, const nvis::vec2& x, const Image<T>& image)
{
    i = (int)floor(image.m * x[0]);
    j = (int)floor(image.n * x[1]);
    image.map_coord(i, j);
}

template< typename T >
bool is_minimum(unsigned int i, unsigned int j, const Image<T>& image, int radius)
{
    T refval = image(i, j);
    int mini = (image.periodic[0] ? i - radius : std::max((int)i - radius, 0));
    int maxi = (image.periodic[0] ? i + radius : std::min((int)i + radius, (int)image.m - 1));
    int minj = (image.periodic[1] ? j - radius : std::max((int)j - radius, 0));
    int maxj = (image.periodic[1] ? j + radius : std::min((int)j + radius, (int)image.n - 1));
    for (int _i = mini ; _i <= maxi ; ++_i) {
        for (int _j = minj ; _j <= maxj ; ++_j) {
            if (_i == i && _j == j) {
                continue;
            }
            T val = image(_i, _j);
            if (val < refval) {
                return false;
            }
        }
    }
    return true;
}

template< typename T >
bool is_maximum(unsigned int i, unsigned int j, const Image<T>& image, int radius)
{
    T refval = image(i, j);
    int mini = (image.periodic[0] ? i - radius : std::max((int)i - radius, 0));
    int maxi = (image.periodic[0] ? i + radius : std::min((int)i + radius, (int)image.m - 1));
    int minj = (image.periodic[1] ? j - radius : std::max((int)j - radius, 0));
    int maxj = (image.periodic[1] ? j + radius : std::min((int)j + radius, (int)image.n - 1));
    for (int _i = mini ; _i <= maxi ; ++_i) {
        for (int _j = minj ; _j <= maxj ; ++_j) {
            if (_i == i && _j == j) {
                continue;
            }
            T val = image(_i, _j);
            if (val > refval) {
                return false;
            }
        }
    }
    return true;
}

// struct Extremum {
//  int i, j;
//  float val;
//  bool max;
// };
//
// struct Lt_Extremum {
//  bool operator()(const Extremum& a, const Extremum& b) const {
//      return a.val < b.val;
//  }
// };

struct Extremum {
    nvis::vec2 pos;
    float val;
    bool max;
};

struct Lt_Extremum {
    bool operator()(const Extremum& a, const Extremum& b) const {
        return a.val < b.val;
    }
};

#endif
























