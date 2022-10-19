#include <iostream>
#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>

double K;
const double twopi = 6.28318530717958647688;

nvis::vec2 next(const nvis::vec2& x)
{
    nvis::vec2 y(x);
    y[0] += K/twopi * sin(twopi*y[1]);
    y[1] += y[0];
}

double modulo(double d)
{
    double m;
    if (d<0) {
        m = -fmod(fabs(d), 1);
    } else {
        m = fmod(d, 1);
    }
    if (m<0) {
        m += 1;
    }
    return m;
}

int main(int argc, char* argv[])
{
    K = atof(argv[1]);
    
    float* data = (float*)calloc(3*1000*1000, sizeof(float));
    
    srand48(time(0));
    for (int j=0 ; j<1000 ; ++j) {
        for (int i=0 ; i<1000 ; ++i) {
            nvis::vec2 x((double)i/1000, (double)j/1000);
            nvis::vec2 y(x);
            nvis::fvec3 col(drand48(), drand48(), drand48());
            for (int n=0 ; n<200 ; ++n) {
                int u = floor(1000*modulo(y[0]));
                int v = floor(1000*modulo(y[1]));
                data[3*(1000*u+v)  ] = col[0];
                data[3*(1000*u+v)+1] = col[1];
                data[3*(1000*u+v)+2] = col[2];
                y = next(y);
            }
        }
    }
    
    Nrrd* nout = nrrdNew();
    size_t dim[3] = {3, 1000, 1000};
    nrrdWrap_nva(nout, data, nrrdTypeFloat, 3, dim);
    nrrdSave(argv[2], nout, NULL);
    
    return 0;
}
