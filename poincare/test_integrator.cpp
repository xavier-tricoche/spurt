#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <tokamak/poincare_map.hpp>
#include <math/dopri5.hpp>
#include <__poincare/__map.hpp>
#include <__poincare/__metric.hpp>
#include <iostream>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <sstream>


nvis::vec3 integrate(const nvis::vec3& init, double delta_t,
                     const tokamak_nimrod_parametric& field)
{
    typedef nvis::dopri5<nvis::vec3>    int_type;
    typedef int_type::step              int_step;
    
    nvis::dopri5<nvis::vec3>            intg;
    
    intg.t = 0;
    intg.y = init;
    intg.t_max = delta_t;
    intg.h = 1e-1;
    intg.abstol = 1.0e-7;
    intg.reltol = 1.0e-7;
    nvis::dopri5<nvis::vec3>::result res;
    nvis::dopri5<nvis::vec3>::step   step;
    
    for (;;) {
        res = intg.do_step(field, step);
        
        // std::cout << "currently at " << step.y1() << "\n";
        
        if (res == int_type::T_MAX_REACHED) {
            break;
        } else if (res != int_type::OK) {
            throw;
        }
    }
    
    return step.y1();
}

nvis::mat3 approximate_jacobian(const nvis::vec3& init, double delta_t,
                                const tokamak_nimrod_parametric& field)
{
    double h = 0.05;
    nvis::vec3 ddx = 0.5 / h * (integrate(init + nvis::vec3(h, 0, 0), delta_t, field) -
                                integrate(init - nvis::vec3(h, 0, 0), delta_t, field));
    nvis::vec3 ddy = 0.5 / h * (integrate(init + nvis::vec3(0, h, 0), delta_t, field) -
                                integrate(init - nvis::vec3(0, h, 0), delta_t, field));
    nvis::vec3 ddz = 0.5 / h * (integrate(init + nvis::vec3(0, 0, h), delta_t, field) -
                                integrate(init - nvis::vec3(0, 0, h), delta_t, field));
                                
    nvis::mat3 J;
    for (int i = 0 ; i < 3 ; ++i) {
        J[i][0] = ddx[i];
        J[i][1] = ddy[i];
        J[i][2] = ddz[i];
    }
    
    return J;
}

nvis::mat3 integrate_jacobian(const nvis::vec3& init, double delta_t,
                              const tokamak_nimrod_parametric& field)
{
    typedef nvis::dopri5<nvis::vec3>    int_type;
    typedef int_type::step              int_step;
    
    nvis::dopri5<nvis::vec3>            intg;
    
    intg.t = 0;
    intg.y = init;
    intg.t_max = delta_t;
    intg.h = 1e-1;
    intg.abstol = 1.0e-7;
    intg.reltol = 1.0e-7;
    nvis::dopri5<nvis::vec3>::result res;
    nvis::dopri5<nvis::vec3>::step   step;
    
    nvis::mat3 J = nvis::mat3::identity();
    
    for (;;) {
        res = intg.do_step(field, step);
        
        // std::cout << "currently at " << step.y1() << "\n";
        
        if (res != int_type::OK && res != int_type::T_MAX_REACHED) {
            throw;
        }
        
        double t0 = step.t0();
        double t1 = step.t1();
        double dt = (t1 - t0) / 10.;
        nvis::mat3 dvdx;
        for (int i = 0 ; i < 10 ; ++i) {
            double t = t0 + (double)i * dt;
            nvis::vec3 x = step.y(t);
            dvdx = field.derivative(0, x);
            J += dt * (dvdx * J);
        }
        
        if (res == int_type::T_MAX_REACHED) {
            break;
        }
    }
    
    return J;
}

int main(int argc, char* argv[])
{
    tokamak_nimrod_parametric field("/scratch2/data/NIMROD/CDXU/hdf5/cdxub7n.h5",
                                    "2000");
                                    
    poincare_map pmap(&field);
    pmap.precision(1.0e-7);
    
    map_metric::bbox = field.bounds();
    map_metric::periodic[0] = true;
    map_metric::periodic[1] = false;
    
    srand48(time(0));
    double error1 = 0., error2 = 0.;
    double maxerr1 = 0., maxerr2 = 0.;
    nvis::vec3 f, f1, f2, df, g;
    nvis::mat3 J1, J2;
    
    int count = 0;
    double avg_error = 0.;
    std::vector< nvis::vec3 >   hits;
    std::vector< nvis::mat3 >   jacs;
    std::vector< double >       times;
    for (unsigned i = 0 ; i < 10 ; ++i) {
        nvis::vec2 x(120.*drand48(), 80.*drand48());
        nvis::vec3 x3d(x[1], x[0], 0);
        
        std::cout << "\nintegrating from " << x << "\n\n\n";
        
        pmap.debug_jmap(x, hits, jacs, times, 20);
        for (unsigned int j = 0 ; j < times.size() ; ++j) {
        
            spurt::map_wrapper< poincare_map > map(pmap, j + 1);
            spurt::central_diff_jacobian< spurt::map_wrapper< poincare_map > > cdj(map, 0.05, 0.05);
            spurt::integral_jacobian< spurt::map_wrapper< poincare_map > > itj(map);
            nvis::mat2 __itj, __cdj;
            try {
                __cdj = cdj(x);
                __itj = itj(x);
            } catch (...) {
                std::cerr << "exception caught at " << x << '\n';
                break;
            }
            
            nvis::mat3 Japp = approximate_jacobian(x3d, times[j], field);
            std::cout << "after " << j + 1 << " rotations:\nintegrated jacobian = " << jacs[j]
                      << "corresponding 2d integrated jacobian = " << __itj
                      << "central difference jacobian = " << Japp
                      << "corresponding 2d central difference jacobian = " << __cdj;
            double err = nvis::norm(jacs[j] - Japp) / nvis::norm(Japp);
            std::cout << "error in 3d = " << err*100 << "%\n\n";
            err = nvis::norm(__itj - __cdj) / nvis::norm(__cdj);
            std::cout << "error in 2d = " << err*100 << "%\n\n";
            avg_error += err;
            ++count;
        }
    }
    avg_error /= (double)count;
    
    std::cout << "average error = " << avg_error*100 << "%\n";
    
    return 0;
}





































































