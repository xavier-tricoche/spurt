#include <iostream>

#include <vector>
#include <string>
#include <algorithm>

#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>

template< typename T >
struct Image {
    Image(Nrrd* nrrd)
            : m(nrrd->axis[0].size), n(nrrd->axis[1].size), data((T*)nrrd->data) {}

    const T& operator()(unsigned int i, unsigned int j) const {
        return data[i+j*m];
    }
    T& operator()(unsigned int i, unsigned int j) {
        return data[i+j*m];
    }
    nvis::vec2 gradient(unsigned int i, unsigned int j) const {
        const Image& self = *this;
        return 0.5*nvis::vec2(self(i + 1, j) - self(i - 1, j),
                              self(i, j + 1) - self(i, j - 1));
    }
    nvis::vec3 hessian(unsigned int i, unsigned int j) const {
        nvis::vec2 dgdx = 0.5 * (gradient(i + 1, j) - gradient(i - 1, j));
        nvis::vec2 dgdy = 0.5 * (gradient(i, j + 1) - gradient(i, j - 1));
        return nvis::vec3(dgdx[0], dgdx[1], dgdy[1]);
    }

    T* data;
    unsigned int m, n;
};


void eigen(double& lmin, nvis::vec2& evmin, const nvis::vec3& H)
{
    double tr = H[0] + H[2];
    double det = H[0] * H[2] - H[1] * H[1];
    /* l^2 - tr*l + det = 0 */
    /* delta = tr*tr - 4*det */
    /* l = 0.5*(tr - sqrt(delta)) */
    double delta = tr * tr - 4.0 * det;
    lmin = 0.5 * (tr - sqrt(delta));
    /*
        (a-l)x +         by = 0
            bx +     (c-l)y = 0
    */
    double a = fabs(H[0] - lmin);
    double b = fabs(H[2] - lmin);
    if (a > b) {
        evmin[0] = -H[1];
        evmin[1] = H[0] - lmin;
    }
    else {
        evmin[0] = H[2] - lmin;
        evmin[1] = -H[1];
    }
    evmin /= nvis::norm(evmin);
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "USAGE: " << argv[0] << " <input.nrrd> <output.nrrd> <threshold>\n";
        exit(-1);
    }
    float threshold = 0;
    if (argc == 4) {
        threshold = atof(argv[3]);
    }

    Nrrd *nin = spurt::readNrrd(argv[1]);
    Image< float > image(nin);
    unsigned int M = image.m;
    unsigned int N = image.n;

    float *data = (float *)nin->data;
    float minval = *std::min_element(&data[0], &data[M*N]);
    float maxval = *std::max_element(&data[0], &data[M*N]);

    float * rgb = (float*)calloc(3 * M * N, sizeof(float));

    std::vector< double > evals;


#pragma openmp parallel for
    for (int n = 0 ; n < image.m*image.n ; ++n) {
        unsigned int i = n % image.m;
        unsigned int j = n / image.n;

        rgb[3*n] = rgb[3*n+1] = rgb[3*n+2] = data[n];

        if (i == 0 || j == 0 || i == image.m - 1 || j == image.n - 1) continue;

        nvis::vec2 emin;
        double lmin;
        nvis::vec3 hess = image.hessian(i, j);
        eigen(lmin, emin, hess);
        if (lmin >= threshold) continue;

        nvis::vec2 grad = image.gradient(i, j);
        double dot = nvis::inner(grad, emin);

        bool found = false;
        for (int di = -1 ; di <= 1 && !found; ++di) {
            for (int dj = -1 ; dj <= 1 && !found; ++dj) {
                if (di == 0 && dj == 0) continue;

                int ii = i + di;
                int jj = j + dj;
                nvis::vec2 g = image.gradient(ii, jj);
                nvis::vec3 H = image.hessian(ii, jj);
                nvis::vec2 e;
                double l;
                eigen(l, e, H);
                if (nvis::inner(e, emin) < 0) {
                    e *= -1;
                }
                double d = nvis::inner(e, g);
                if (d*dot <= 0 && lmin<l) {
                    rgb[3*n] = maxval;
                    rgb[3*n+1] = rgb[3*n+2] = minval;
                    found = true;
                    evals.push_back(lmin);
                }
            }
        }
    }

    std::vector< size_t > size(3);
    size[0] = 3;
    size[1] = M;
    size[2] = N;
    std::vector< float > spacing;
    spurt::writeNrrd((void*)rgb, argv[2], nrrdTypeFloat, size, spacing);

    std::sort(evals.begin(), evals.end());
    std::cout << "ridge strength ranges from " << evals[0] << " and " << evals.back() << '\n';

    return 0;
}





































































