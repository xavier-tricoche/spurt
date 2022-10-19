#include <image/creaseline2d.hpp>
#include <vector>
#include <iostream>
#include <data/raster.hpp>
#include <teem/nrrd.h>
#include <map>
#include <misc/sort.hpp>

Nrrd* readNrrd(const std::string& filename)
{
    Nrrd* nin = nrrdNew();
    
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        char* err = biffGetDone(NRRD);
        std::cerr << "FDTICoherence: " << err << std::endl;
        exit(-1);
    }
    
    return nin;
}

typedef std::pair< unsigned int, double > value;
struct SortValue {
    bool operator()(const value& v0, const value& v1) const {
        return v0.second < v1.second;
    }
};

int main(int argc, char* argv[])
{

    if (argc != 4) {
        std::cout << "USAGE: " << argv[0] << " <in> <distances> <out>\n";
        return -1;
    }
    
    std::string filename = argv[1];
    Nrrd* nin = readNrrd(filename);
    
    Nrrd* _dist = readNrrd(argv[2]);
    float* dist = (float*)_dist->data;
    
    xavier::Raster::Grid grid;
    grid.nx = nin->axis[0].size;
    grid.ny = nin->axis[1].size;
    grid.minx = grid.miny = 0;
    grid.maxx = nin->axis[0].size - 1;
    grid.maxy = nin->axis[1].size - 1;
    
    std::vector< nvis::vec2 > points;
    std::vector< std::list< unsigned int > > ridges;
    float* input = (float*)nin->data;
    unsigned int N = nin->axis[0].size * nin->axis[1].size;
    std::vector< float > tmp(&input[0], &input[N]);
    float max = *std::max_element(tmp.begin(), tmp.end());
    float min = *std::min_element(tmp.begin(), tmp.end());
    float avg = 0;
    for (unsigned int n = 0 ; n < N ; ++n) {
        avg += input[n];
    }
    avg /= (float)N;
    
    Nrrd* _nin = nrrdNew();
    float* _data = (float*)calloc(N * 2, sizeof(float));
    for (unsigned int i = 0 ; i < N ; ++i) {
        _data[i] = _data[i+N] = input[i];
    }
    nrrdWrap_va(_nin, _data, nrrdTypeFloat, 3, nin->axis[0].size, nin->axis[1].size, 2);
    nrrdAxisInfoSet_va(_nin, nrrdAxisInfoCenter, nrrdCenterCell, nrrdCenterCell, nrrdCenterCell);
    nrrdAxisInfoSet_va(_nin, nrrdAxisInfoSpacing, 1.0, 1.0, 1.0);
    xavier::crease::extract(_nin, grid, -1000000, true, points, ridges);
    
    unsigned int pmax = _dist->axis[2].size;
    float* data = (float*)calloc(3 * pmax * N, sizeof(float));
    for (unsigned int i = 0 ; i < N ; ++i) {
        for (unsigned int p = 0 ; p < pmax ; ++p) {
            data[3*N*p+3*i  ] = input[i];
            data[3*N*p+3*i+1] = input[i];
            data[3*N*p+3*i+2] = input[i];
        }
    }
    
    srand48(time(0));
    std::vector< double > ridge_strength(ridges.size(), 0);
    for (unsigned n = 0 ; n < ridges.size() ; ++n) {
        std::list< unsigned int >::iterator it;
        if (ridges[n].size() < 2) {
            // std::cout << "ridge length is only " << ridges[n].size() << "\n";
            continue;
        }
        
        double str = 0;
        for (it = ridges[n].begin() ; it != ridges[n].end() ; ++it) {
            str += fabs(xavier::crease::crease_strength[*it]);
        }
        ridge_strength[n] = str / (double)ridges[n].size();
    }
    std::vector< double > _tmp(ridge_strength.begin(), ridge_strength.end());
    std::vector< unsigned int > sorted;
    xavier::sort(_tmp, sorted);
    double threshold = _tmp[sorted[_tmp.size()/2]];
    
    for (unsigned int p = 0 ; p < pmax ; ++p) {
    
        std::vector< value > ridge_min;
        for (unsigned n = 0 ; n < ridges.size() ; ++n) {
            if (ridge_strength[n] < threshold || ridges[n].size() < 10) {
                continue;
            }
            
            nvis::fvec3 color(0.3 * drand48(), 0.3 * drand48(), 0.3 * drand48());
            color *= max;
            std::list< unsigned int >::iterator it;
            value min = value(0, std::numeric_limits<double>::max());
            for (it = ridges[n].begin() ; it != ridges[n].end() ; ++it) {
            
                const nvis::vec2& x = points[*it];
                unsigned int i = (unsigned int)floor(x[0]);
                unsigned int j = (unsigned int)floor(x[1]);
                unsigned int q = i + j * nin->axis[0].size;
                data[3*N*p+3*q] = color[0];
                data[3*N*p+3*q+1] = color[1];
                data[3*N*p+3*q+2] = color[2];
                
                double d = dist[p*N+q];
                if (d < min.second) {
                    min = value(q, d);
                }
            }
            ridge_min.push_back(min);
        }
        std::sort(ridge_min.begin(), ridge_min.end(), SortValue());
        
        for (unsigned int k = 0 ; k < 2*(p + 1) ; ++k) {
            unsigned int q = ridge_min[k].first;
            unsigned int shift = 3 * N * p;
            unsigned int col = _dist->axis[0].size;
            
            int qs[] = { q - col - 1, q - col, q - col + 1, q - 1, q, q + 1, q + col - 1, q + col, q + col + 1};
            for (unsigned int n = 0 ; n < 9 ; ++n) {
                int x = qs[n];
                if (x < 0) {
                    continue;
                }
                unsigned int i = x % col;
                unsigned int j = x / col;
                if (i >= col || j >= col) {
                    continue;
                }
                
                data[shift+3*x] = 5 * max;
                data[shift+3*x+1] = 0;
                data[shift+3*x+2] = 0;
            }
        }
        
    }
    
    Nrrd* nout = nrrdNew();
    if (nrrdWrap_va(nout, data, nrrdTypeFloat, 4, 3, nin->axis[0].size, nin->axis[1].size, pmax) ||
            nrrdSave(argv[3], nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                  << std::endl;
    }
    nrrdNuke(nout);
    
    nrrdNuke(nin);
    
    return 0;
}



























































