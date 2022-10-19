#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

struct data {
    size_t n;
    double * c;
    double * theta;
    double * sigma;
};

int model_f(const gsl_vector * x, void *data,
            gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *c = ((struct data *)data)->c;
    double *theta = ((struct data *)data)->theta;
    double *sigma = ((struct data *)data)->sigma;

    double iso = gsl_vector_get (x, 0);
    double dA = gsl_vector_get (x, 1);
    double phi = gsl_vector_get (x, 2);

    for (size_t i = 0; i < n; i++)
    {
        double Ci = iso + dA * cos( 2*theta[i] + phi);
        gsl_vector_set (f, i, (Ci - c[i])/sigma[i]);
    }

    return GSL_SUCCESS;
}

int model_df(const gsl_vector * x, void *data,
             gsl_matrix * J)
{
    size_t n = ((struct data *)data)->n;
    double *c = ((struct data *)data)->c;
    double *theta = ((struct data *)data)->theta;
    double *sigma = ((struct data *)data)->sigma;

    double iso = gsl_vector_get (x, 0);
    double dA = gsl_vector_get (x, 1);
    double phi = gsl_vector_get (x, 2);

    for (size_t i = 0; i < n; i++)
    {
        double dCi_diso = 1.0;
        double dCi_ddA = cos( 2*theta[i] + phi);
        double dCi_dphi = -dA * sin (2*theta[i]+phi);
        double s = sigma[i];
        gsl_matrix_set (J, i, 0, dCi_diso/s);
        gsl_matrix_set (J, i, 1, dCi_ddA/s);
        gsl_matrix_set (J, i, 2, dCi_dphi/s);
    }
    return GSL_SUCCESS;
}

int model_fdf (const gsl_vector * x, void *data,
               gsl_vector * f, gsl_matrix * J)
{
    model_f (x, data, f);
    model_df (x, data, J);

    return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
    printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
            "|f(x)| = %g\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_blas_dnrm2 (s->f));
}

int main(int argc, char* argv[])
{
    std::string filename = argv[1];
    std::fstream in(filename.c_str(), std::ios::in); 
    std::vector<double> thetas;
    std::vector<double> cs;
    // skip header
    std::string header;
    in >> header;
    while (!in.eof()) {
        double t, c=-1;
        char s;
        in >> t >> s >> c;
        if (c < 0) break;
        thetas.push_back(fmod(t, M_PI));
        cs.push_back(c);
    }
    in.close();
    size_t N = cs.size();
    
    // initialization
    double iso_guess = 0;
    for (int i=0 ; i<N ; ++i) {
        iso_guess += cs[i];
    }
    iso_guess /= (double)N;
    int max_id = std::distance(cs.begin(), std::max_element(cs.begin(), cs.end()));
    double dA_guess = cs[max_id] - iso_guess;
    double phi_guess = -2*thetas[max_id];
    
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    unsigned int i, iter = 0;
    size_t n = N;
    const size_t p = 3;

    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double c[N], theta[N], sigma[N];
    struct data d = { n, c, theta, sigma};
    gsl_multifit_function_fdf f;
    double x_init[3] = { iso_guess, dA_guess, phi_guess }; //<= for Xavier: initialization here
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    const gsl_rng_type * type;
    gsl_rng * r;

    gsl_rng_env_setup();

    type = gsl_rng_default;
    r = gsl_rng_alloc (type);

    f.f = &model_f;
    f.df = &model_df;
    f.fdf = &model_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    /* This is the data to be fitted */

    for (i = 0; i < n; i++)
    {    
        // for Xavier: Fill data here
        c[i] = cs[i];
        theta[i] = thetas[i];
        sigma[i] = 0.1; // may be smaller or larger accuracy?

        printf ("data: %g %g %g\n", c[i], theta[i], sigma[i]);
    };

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    print_state (iter, s);

    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);

        printf ("status = %s\n", gsl_strerror (status));

        print_state (iter, s);

        if (status)
            break;

        status = gsl_multifit_test_delta (s->dx, s->x,
                                          1e-4, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 500);

    gsl_multifit_covar (s->J, 0.0, covar);
    {
        double chi = gsl_blas_dnrm2(s->f);
        double dof = n - p;
        double k = GSL_MAX_DBL(1, chi / sqrt(dof));

        printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

        printf ("iso      = %.5f +/- %.5f\n", FIT(0), k*ERR(0));
        printf ("dA       = %.5f +/- %.5f\n", FIT(1), k*ERR(1));
        printf ("phi      = %.5f +/- %.5f\n", FIT(2), k*ERR(2));
    }

    printf ("status = %s\n", gsl_strerror (status));

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
    return 0;
}