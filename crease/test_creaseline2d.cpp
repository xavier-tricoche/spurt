#include <image/creaseline2d.hpp>
#include <vector>
#include <iostream>
#include <data/raster.hpp>
#include <teem/nrrd.h>

Nrrd* readNrrd(const std::string& filename)
{
    Nrrd *nin = nrrdNew();

    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        char *err = biffGetDone(NRRD);
        std::cerr << "FDTICoherence: " << err << std::endl;
        exit(-1);
    }

    return nin;
}

int main(int argc, char* argv[])
{

    if (argc != 3) {
        std::cout << "USAGE: " << argv[0] << " <in> <out>\n";
        return -1;
    }

    std::string filename = argv[1];
    Nrrd *nin = readNrrd(filename);

    spurt::Raster::Grid grid;
    grid.nx = nin->axis[0].size;
    grid.ny = nin->axis[1].size;
    grid.minx = grid.miny = 0;
    grid.maxx = nin->axis[0].size - 1;
    grid.maxy = nin->axis[1].size - 1;

    std::vector< nvis::vec2 > points;
    std::vector< std::list< unsigned int > > ridges;
    float *input = (float*)nin->data;
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
    float *_data = (float*)calloc(N * 2, sizeof(float));
    for (unsigned int i = 0 ; i < N ; ++i) {
        _data[i] = _data[i+N] = input[i];
    }
    nrrdWrap_va(_nin, _data, nrrdTypeFloat, 3, nin->axis[0].size, nin->axis[1].size, 2);
    nrrdAxisInfoSet_va(_nin, nrrdAxisInfoCenter, nrrdCenterCell, nrrdCenterCell, nrrdCenterCell );
    nrrdAxisInfoSet_va(_nin, nrrdAxisInfoSpacing, 1.0, 1.0, 1.0 );
    spurt::crease::extract(_nin, grid, -1000000, true, points, ridges);
    float *data = (float*)calloc(3 * N, sizeof(float));
    for (unsigned int i = 0 ; i < N ; ++i) {
        data[3*i  ] = input[i];
        data[3*i+1] = input[i];
        data[3*i+2] = input[i];
    }

    /*
    for (unsigned int n = 0 ; n < spurt::crease::skipped.size() ; ++n) {
        if (spurt::crease::skipped[n]) {
            unsigned int i = n % grid.nx;
            unsigned int j = n / grid.nx;

            unsigned int ii = i * nin->axis[0].size / grid.nx;
            unsigned int jj = j * nin->axis[1].size / grid.ny;

            unsigned int k = ii + jj * nin->axis[0].size;
            data[3*k  ] = 0;
            data[3*k+1] = 0;
            data[3*k+2] = max;
        }
    }
    */
    /*
    for (unsigned int n = 0 ; n < spurt::crease::singular.size() ; ++n) {
        if (spurt::crease::singular[n]) {
            unsigned int i = n % grid.nx;
            unsigned int j = n / grid.nx;

            unsigned int ii = i * nin->axis[0].size / grid.nx;
            unsigned int jj = j * nin->axis[1].size / grid.ny;

            unsigned int k = ii + jj * nin->axis[0].size;
            data[3*k  ] = 0;
            data[3*k+1] = max;
            data[3*k+2] = 0;
        }
    }
    */

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
            str += fabs(spurt::crease::crease_strength[*it]);
        }
        ridge_strength[n] = str / (double)ridges[n].size();
    }
    std::vector< double > _tmp(ridge_strength.begin(), ridge_strength.end());
    std::sort(_tmp.begin(), _tmp.end());
    double threshold = _tmp[_tmp.size()/2];

    for (unsigned n = 0 ; n < ridges.size() ; ++n) {
        if (ridge_strength[n] < threshold || ridges[n].size()<10 ) continue;

        nvis::fvec3 color(0.3 * drand48(), 0.3 * drand48(), 0.3 * drand48());
        color *= max;
        std::list< unsigned int >::iterator it;
        for (it = ridges[n].begin() ; it != ridges[n].end() ; ++it) {

            const nvis::vec2& x = points[*it];
            unsigned int i = (unsigned int)floor(x[0]);
            unsigned int j = (unsigned int)floor(x[1]);
            unsigned int p = i + j * nin->axis[0].size;
            data[3*p] = color[0];
            data[3*p+1] = color[1];
            data[3*p+2] = color[2];
        }
    }

    Nrrd *nout = nrrdNew();
    if (nrrdWrap_va(nout, data, nrrdTypeFloat, 3, 3, nin->axis[0].size, nin->axis[1].size) ||
        nrrdSave(argv[2], nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
        << std::endl;
    }
    nrrdNuke(nout);

    nrrdNuke(nin);

    return 0;
}




































