#ifndef __CREASE_HPP__
#define __CREASE_HPP__

#include <math/fixed_vector.hpp>
#include <image/probe.hpp>
#include <image/eigen.hpp>
#include <vector>
#include <list>

namespace xavier {
namespace crease {
// this parameter controls the subdivision going on in the
// detection of crease points on provided edges.
// 0 value will lead to no subdivision
unsigned int nsub;

bool find_intersection(double* x0, double* x1, double* inter,
                       gage_interface::scalar_wrapper& gH_wrapper,
                       bool ridge);
bool find_intersection(const nvis::vec2& x0, const nvis::vec2& x1,
                       nvis::vec2& inter,
                       gage_interface::scalar_wrapper& gH_wrapper,
                       bool ridge);
                       
bool find_intersection(const nvis::vec2& x0, const nvis::vec2& x1,
                       nvis::vec2& inter,
                       gage_interface::scalar_wrapper& gH_wrapper,
                       bool ridge, double& strength);
                       
bool quick_search(double* x0, double* x1, double* inter,
                  gage_interface::scalar_wrapper& gH,
                  bool ridge);
bool quick_search(const nvis::vec2& x0, const nvis::vec2& x1,
                  nvis::vec2& inter,
                  gage_interface::scalar_wrapper& gH,
                  bool ridge);
                  
bool eigen(double* emin, double& lmin, double* H,
           bool ridge);
bool eigen(nvis::vec2& emin, double& lmin, const nvis::vec3& H,
           bool ridge);
};
};

bool xavier::crease::
find_intersection(double* x0, double* x1, double* inter,
                  gage_interface::scalar_wrapper& gH_wrapper,
                  bool ridge)
{
    if (nsub == 0)
        return quick_search(x0, x1, inter, gH_wrapper, ridge);

    // subdivide edge
    double x[nsub+2][2];
    for (unsigned int i = 0 ; i < nsub + 2 ; i++) {
        double a = (double)i / (double)(nsub + 1);
        x[i][0] = (1. - a) * x0[0] + a * x1[0];
        x[i][1] = (1. - a) * x0[1] + a * x1[1];
    }

    // compute gradient and reference eigenvector at one vertex
    double evec[nsub+2][2];
    double dot0, dot1;
    double H[3], g[2];
    double lambda;
    std::vector< double > roots;
    gH_wrapper.hessian(x[0], H);
    bool valid0 = eigen(evec[0], lambda, H, ridge);
    if (valid0) {
        if (!gH_wrapper.gradient(x[0], g)) return false;
        dot0 = evec[0][0] * g[0] + evec[0][1] * g[1];
    }
    else {
        std::cout << "false #1\n";
        return false;
    }
    for (unsigned int i = 1 ; i < nsub + 2 ; i++) {
        if (!gH_wrapper.hessian(x[i], H)) return false;
        bool valid1 = eigen(evec[i], lambda, H, ridge);
        double _dot = evec[i-1][0] * evec[i][0] + evec[i-1][1] * evec[i][1];
        if (_dot < 0) {
            evec[i][0] *= -1;
            evec[i][1] *= -1;
            // would need a recursive subdivision if dot product is
            // lower than a prescribed value that ensures the strong
            // similarity between consecutive eigenvector directions
            // but oh well...
        }

        if (valid1) {
            if (!gH_wrapper.gradient(x[i], g)) {
                std::cout << "false #2\n";
                return false;
            }
            dot1 = evec[i][0] * g[0] + evec[i][1] * g[1];
            if (valid0 && dot1*dot0 < 0) {
                double t = -dot0 / (dot1 - dot0);
                roots.push_back((1. - t)*x[i-1][0] + t*x[i][0]);
                roots.push_back((1. - t)*x[i-1][1] + t*x[i][1]);
            }

            dot0 = dot1;
            valid0 = true;
        }
        else
        {
            std::cout << "false #3\n";
            return false;
        }
    }

    unsigned int nbroots = roots.size() / 2;

    std::cout << nbroots << " roots found" << std::endl;

    if (!nbroots) return false;
    else if (nbroots == 1) {
        inter[0] = roots[0];
        inter[1] = roots[1];
        return true;
    }
    else if (nbroots % 2) {
        inter[0] = inter[1] = 0;
        for (unsigned int n = 0 ; n < nbroots ; n++) {
            inter[0] += roots[2*n  ];
            inter[1] += roots[2*n+1];
        }
        inter[0] /= (double)nbroots;
        inter[1] /= (double)nbroots;
        return true;
    }
    else return false;
}

bool xavier::crease::
quick_search(double* x0, double* x1, double* inter,
             gage_interface::scalar_wrapper& gH,
             bool ridge)
{
    // compute gradient and reference eigenvector at one vertex
    double evec[2][2];
    double dot0, dot1;
    double H[3], g[2];
    double lambda;
    // get relevant eigenvectors and check eigenvalues
    if (!gH.hessian(x0, H) ||
        !eigen(evec[0], lambda, H, ridge)) return false;
    if (!gH.hessian(x1, H) ||
        !eigen(evec[1], lambda, H, ridge)) return false;

    // enforce orientation consistency across eigenvectors
    if (evec[0][0]*evec[1][0] + evec[0][1]*evec[1][1] < 0) {
        evec[1][0] *= -1;
        evec[1][1] *= -1;
    }

    // check for zero crossing of dot product with gradient
    if (!gH.gradient(x0, g)) return false;
    dot0 = evec[0][0] * g[0] + evec[0][1] * g[1];
    if (!gH.gradient(x1, g)) return false;
    dot1 = evec[1][0] * g[0] + evec[1][1] * g[1];

    if (dot0*dot1 <= 0) {
        if (dot0 == dot1) {
            inter[0] = x0[0];
            inter[1] = x0[1];
        }
        else {
            double u = -dot0 / (dot1 - dot0);
            inter[0] = (1 - u) * x0[0] + u * x1[0];
            inter[1] = (1 - u) * x0[1] + u * x1[1];
        }

        return true;
    }

    return false;
}

bool xavier::crease::
eigen(double* evec, double& lambda, double* H, bool ridge)
{
    static nvis::vec3 Hessian;
    static nvis::vec2 emin, emax;
    double lmin, lmax;

    Hessian[0] = H[0];
    Hessian[1] = H[1];
    Hessian[2] = H[2];

    xavier::eigensystem(emin, emax, lmin, lmax, Hessian);
    if (ridge) {
        lambda = lmin;
        evec[0] = emin[0];
        evec[1] = emin[1];
        return lmin < 0;
    }
    else {
        lambda = lmax;
        evec[0] = emax[0];
        evec[1] = emax[1];
        return lmax > 0;
    }
}

bool xavier::crease::
find_intersection(const nvis::vec2& x0, const nvis::vec2& x1,
                  nvis::vec2& inter,
                  gage_interface::scalar_wrapper& gH_wrapper,
                  bool ridge, double& strength)
{
    strength = 0;
    if (find_intersection(x0, x1, inter, gH_wrapper, ridge)) {
        nvis::vec2 evec;
        double eval;
        nvis::vec3 H;
        nvis::vec3 p(inter[0], inter[1], 0.5);
        gH_wrapper.hessian(inter, H);

//       std::cout << "Hessian at " << inter << " is " << H << std::endl;

        if (eigen(evec, eval, H, ridge)) {
//       std::cout << "checking lambda (" << eval << ")" << std::endl;
//       nvis::vec2 He( H[0]*evec[0]+H[1]*evec[1],
//              H[1]*evec[0]+H[2]*evec[1] );
//       std::cout << "cross product: " << fabs( nvis::cross( He, evec ) )
//             << ", lambda test: " << nvis::norm( He-eval*evec )
//             << std::endl;

//       nvis::vec3 evals;
//       gH_wrapper.hess_evals( p, evals );
//       if ( ridge )
//         std::cout << "min Hess_Eval = " << evals[2] << std::endl;
//       else
//         std::cout << "max Hess_Eval = " << evals[0] << std::endl;


            strength = eval;
            return true;
        }
    }

    return false;
}

bool xavier::crease::
find_intersection(const nvis::vec2& x0, const nvis::vec2& x1,
                  nvis::vec2& inter,
                  gage_interface::scalar_wrapper& gH_wrapper,
                  bool ridge)
{
    if (nsub == 0)
        return quick_search(x0, x1, inter, gH_wrapper, ridge);

    // subdivide edge
    std::vector< nvis::vec2 > x(nsub + 2);
    for (unsigned int i = 0 ; i < nsub + 2 ; i++) {
        double a = (double)i / (double)(nsub + 1);
        x[i] = (1. - a) * x0 + a * x1;
    }

    // check Hessian's eigenvalues at end vertices
    std::vector< nvis::vec2 > evec(nsub + 2);
    nvis::vec3 H;
    nvis::vec2 g;
    double lambda;
    std::vector< nvis::vec2 > roots;
    if (!gH_wrapper.hessian(x[0], H) ||
        !eigen(evec[0], lambda, H, ridge))
        return false;
    if (!gH_wrapper.hessian(x[nsub+1], H) ||
        !eigen(evec[nsub+1], lambda, H, ridge))
        return false;

    // loop over segment discretization
    for (unsigned int i = 0 ; i < nsub + 1 ; i++) {
        if (quick_search(x[i], x[i+1], inter, gH_wrapper, ridge))
            return true;
    }

    return false;
}

bool xavier::crease::
quick_search(const nvis::vec2& x0, const nvis::vec2& x1, nvis::vec2& inter,
             gage_interface::scalar_wrapper& gH,
             bool ridge)
{
    // compute gradient and reference eigenvector at one vertex
    std::vector< nvis::vec2 > evec(2);
    double dot0, dot1;
    nvis::vec3 H;
    nvis::vec2 g;
    double lambda;
    // get relevant eigenvectors and check eigenvalues
    if (!gH.hessian(x0, H) ||
        !eigen(evec[0], lambda, H, ridge) ||
        !gH.hessian(x1, H) ||
        !eigen(evec[1], lambda, H, ridge))
        return false;

    // enforce orientation consistency across eigenvectors
    if (nvis::inner(evec[0], evec[1]) < 0) {
        evec[1] *= -1;
    }

    // check for zero crossing of dot product with gradient
    if (!gH.gradient(x0, g)) return false;
    dot0 = nvis::inner(evec[0], g);
    if (!gH.gradient(x1, g)) return false;
    dot1 = nvis::inner(evec[1], g);

    if (dot0*dot1 <= 0) {
        if (dot0 == dot1) {
            inter = x0;
        }
        else {
            double u = -dot0 / (dot1 - dot0);
            inter = (1 - u) * x0 + u * x1;
        }
        return true;
    }

    return false;
}

bool xavier::crease::
eigen(nvis::vec2& evec, double& lambda, const nvis::vec3& H, bool ridge)
{
    nvis::vec2 emin, emax;
    double lmin, lmax;

    xavier::eigensystem(emin, emax, lmin, lmax, H);

    /*
    // checking results
    std::cout << "checking eigen" << std::endl;
    double Hnorm = sqrt( H[0]*H[0] + H[1]*H[1] + H[2]*H[2] );
    double det = ( H[0]-lmin )*( H[2]-lmin )-H[1]*H[1];
    std::cout << "det( H - lmin*I )/||H|| = " << det/Hnorm << std::endl;
    det = ( H[0]-lmax )*( H[2]-lmax )-H[1]*H[1];
    std::cout << "det( H - lmax*I )/||H|| = " << det/Hnorm << std::endl;

    nvis::vec2 Hemin( ( H[0]-lmin )*emin[0]+H[1]*emin[1],
              H[1]*emin[0]+( H[2]-lmin )*emin[1] );
    std::cout << "||( H-lmin*I )emin|| = " << nvis::norm( Hemin ) << std::endl;
    nvis::vec2 Hemax( ( H[0]-lmax )*emax[0]+H[1]*emax[1],
              H[1]*emax[0]+( H[2]-lmax )*emax[1] );
    std::cout << "||( H-lmax*I )emax|| = " << nvis::norm( Hemax ) << std::endl;

    */

    if (ridge) {
        lambda = lmin;
        evec = emin;
        return lmin < 0;
    }
    else {
        lambda = lmax;
        evec = emax;
        return lmax > 0;
    }
}

#endif


