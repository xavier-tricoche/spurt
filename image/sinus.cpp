#include <teem/nrrd.h>
#include <iostream>
#include <vector>
#include <math/fixed_vector.hpp>

const unsigned int size = 5;

void gaussian(int i, int j, int N, float* raster)
{
    for (int jj = j - size ; jj < j + size ; ++jj) {
        for (int ii = i - size ; ii < i + size ; ++ii) {
            if (ii < 0 || ii >= N || jj < 0 || jj >= N) continue;

            float dist = (i - ii) * (i - ii) + (j - jj) * (j - jj);
            unsigned int n = ii + jj * N;
            raster[n] += exp(-dist / 10);
        }
    }

}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "USAGE: " << argv[0] << " <size> <filename>\n";
        return -1;
    }

    int N = atoi(argv[1]);

    float *data = (float *)calloc(2 * N * N, sizeof(float));

    for (unsigned int i = 0 ; i < N ; ++i) {

        float x = - M_PI + (double)i/(double)N * 2.0 * M_PI;
        float y = sin(x);

        unsigned int j = N / 2 + (int)floor(y * N / 3.0f);
        gaussian(i, j, N, data);
    }
    for (unsigned int n = 0 ; n < N ; ++n) {
        data[n+N] = data[n];
    }

    Nrrd *nout = nrrdNew();
    nrrdWrap_va(nout, data, nrrdTypeFloat, 3, N, N, 2);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindSpace, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterCell, nrrdCenterCell, nrrdCenterCell);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, 1.0, 1.0, 1.0);
    if (nrrdSave(argv[2], nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
        << std::endl;
    }

    return 0;

}










