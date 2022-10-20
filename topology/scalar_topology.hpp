#ifndef __SCALAR_TOPOLOGY_HPP__
#define __SCALAR_TOPOLOGY_HPP__

#include <array>
#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include <Eigen/Core>

#include <math/fixed_vector.hpp>
#include <math/angle.hpp>
#include <math/vector_manip.hpp>
#include <image/probe.hpp>
#include <image/nrrd_wrapper.hpp>
#include <image/nrrd_field.hpp>
#include <misc/log_helper.hpp>
#include <misc/progress.hpp>
#include <data/raster.hpp>
#include <vtk/vtk_data_helper.hpp>
#include <vtk/vtk_io_helper.hpp>


namespace {
    template<typename T, int N>
    inline std::string flatten(const Eigen::Matrix<T, N, 1>& v) {
        std::ostringstream oss;
        oss << "[";
        for (int i=0; i<N-1; ++i) {
            oss << v[i] << ", ";
        }
        oss << v[N-1] << "]";
        return oss.str();
    }

} // anonymous


namespace spurt { namespace topology { namespace scalar {

template<typename Matrix_, typename Vector_>
struct naive_solver {
};

template<typename T, int N>
struct naive_solver< Eigen::Matrix<T, N, N >, Eigen::Matrix<T, N, 1> > {
    typedef Eigen::Matrix<T, N, N > matrix_t;
    typedef Eigen::Matrix<T, N, 1> vector_t;
    
    static vector_t solve(const matrix_t& A, const vector_t& b) {
        // std::cout << "solver:\n";
        // std::cout << "        A=\n" << A << "\n";
        // std::cout << "        b=" << flatten(b) << '\n';
        // std::cout << "        inverse(A)=\n" << A.inverse() << '\n';
        // std::cout << "        inverse(A)*A=\n" << A.inverse()*A << '\n';
        // vector_t x = A.inverse()*b;
        // vector_t bb = A*x;
        // std::cout << "        inverse(A)*b=" << flatten(x) << '\n';
        // std::cout << "        Ax=" << flatten(bb) << '\n';
        return A.inverse()*b;
    }
};

template<typename Scalar_, size_t N, 
         typename Vector_ = Eigen::Matrix<Scalar_, N, 1>,
         typename Matrix_ = Eigen::Matrix<Scalar_, N, N>,
         typename Solver_ = naive_solver<Matrix_, Vector_> >
struct SmoothScalarField {
    
    constexpr static size_t dim = N;
    typedef Scalar_ value_t;
    typedef Vector_ pos_t;
    typedef Vector_ gradient_t;
    typedef Matrix_ hessian_t;
    typedef spurt::gage_interface::scalar_wrapper wrapper_t;
    typedef std::pair<pos_t, pos_t> bounds_t;
    typedef Solver_ solver_t;
    
    std::string error_msg(const pos_t& p) const {
        std::ostringstream os;
        os << "Unable to interpolate scalar field " << m_name << " at " << flatten(p) << '\n';
        return os.str();
    }
    
    const bounds_t& bounds() const {
        return m_bounds;
    }

    SmoothScalarField(Nrrd* nin, const std::string& name="anonymous") 
        : m_wrapper(nrrd_utils::make3d<value_t>(nin, true), 
                    spurt::gage_interface::BSPL3_INTERP, 
                    true, true, true, true), m_name(name) {
        m_wrapper.use_world();
        std::vector< std::pair< double, double > > bounds;
        nrrd_utils::compute_raster_bounds(bounds, nin, true);
        for (int i=0; i<dim; ++i) {
            const std::pair<double, double>& p = bounds[i];
            m_bounds.first [i] = p.first;
            m_bounds.second [i] = p.second;
        }
    }

    ~SmoothScalarField() {}
    
    value_t operator()(const pos_t& x) const {
        value_t s;
        if (!m_wrapper.value<pos_t, 2, value_t>(x, s)) {
            throw std::runtime_error(error_msg(x));
        }
        return s;
    }
    
    gradient_t gradient(const pos_t& x) const {
        gradient_t g;
        if (!m_wrapper.gradient<pos_t, 2, gradient_t>(x, g)) {
            throw std::runtime_error(error_msg(x));
        }
        return g;
    }
    
    hessian_t hessian(const pos_t& x) const {
        hessian_t H;
        if (!m_wrapper.hessian<pos_t, 2, hessian_t>(x, H)) {
            throw std::runtime_error(error_msg(x));
        }
        return H;
    }
    
    void hessian_eigen(const pos_t& x, value_t& lmax, value_t& lmin, 
                       gradient_t& emax, gradient_t& emin) const {
        gradient_t evals;
        std::vector< gradient_t > evecs;
        if (!m_wrapper.hess_evals<pos_t, 2, gradient_t>(x, evals)) {
            throw std::runtime_error("Unable to compute Hessian eigenvalues at " + flatten(x));
        }
        std::cout << "computed eigenvalues\n";
        if (!m_wrapper.hess_evecs<pos_t, 2, gradient_t>(x, evecs)) {
            throw std::runtime_error("Unable to compute Hessian eigenvectors at " + flatten(x));
        }  
        std::cout << "computed eigenvectors\n";
        if (evals[0] < evals[1]) {
            lmax = evals[1];
            lmin = evals[0];
            emax = evecs[1];
            emin = evecs[0];
        }   
        else {
            lmax = evals[0];
            lmin = evals[1];
            emax = evecs[0];
            emin = evecs[1];
        }          
    }

    wrapper_t m_wrapper;
    std::string m_name;
    bounds_t m_bounds;
};


template<typename Scalar_, size_t N,
         typename Vector_ = Eigen::Matrix<double, N, 1>,
         typename Matrix_ = Eigen::Matrix<double, N, N>,
         typename Solver_ = naive_solver<double, Vector_> >
struct PWLScalarField {
    constexpr static size_t dim = N;
    typedef Scalar_ value_t;
    typedef Vector_ pos_t;
    typedef Vector_ gradient_t;
    typedef Matrix_ hessian_t;
    typedef spurt::nrrd_field<value_t, N, value_t> field_t;
    typedef std::pair<pos_t, pos_t> bounds_t;
    typedef Solver_ solver_t;
    
    std::string error_msg(const pos_t& p) const {
        std::ostringstream os;
        os << "Unable to interpolate scalar field " << m_name << " at " << flatten(p) << '\n';
        return os.str();
    }
    
    const bounds_t& bounds() const {
        return m_field.grid().bounds();
    }
    
    PWLScalarField(Nrrd* nin, const std::string& name="anonymous") 
    : m_field(nin), m_name(name) {
        m_default_h = 0.05*m_field.rid().spacing();
    }

    ~PWLScalarField() {}
    
    value_t operator()(const pos_t& x) const {
        return m_field.interpolate(x);
    }
    
    gradient_t gradient(const pos_t& x) const {
        gradient_t g;
        auto _g = m_field.differentiate(x);
        for (int i=0; i<N; ++i) g[i] = _g[i];
        return g;
    }
    
    // hessian_t hessian(const pos_t& x) const {
    //     hessian_t H;
    //
    //     return H;
    // }
    //
    // void hessian_eigen(const pos_t& x, value_t& lmax, value_t& lmin,
    //                    gradient_t& emax, gradient_t& emin) const {
    //     gradient_t evals;
    //     std::vector< gradient_t > evecs;
    //     if (!m_wrapper.hess_evals<pos_t, 2, gradient_t>(x, evals)) {
    //         throw std::runtime_error("Unable to compute Hessian eigenvalues at " + flatten(x));
    //     }
    //     std::cout << "computed eigenvalues\n";
    //     if (!m_wrapper.hess_evecs<pos_t, 2, gradient_t>(x, evecs)) {
    //         throw std::runtime_error("Unable to compute Hessian eigenvectors at " + flatten(x));
    //     }
    //     std::cout << "computed eigenvectors\n";
    //     if (evals[0] < evals[1]) {
    //         lmax = evals[1];
    //         lmin = evals[0];
    //         emax = evecs[1];
    //         emin = evecs[0];
    //     }
    //     else {
    //         lmax = evals[0];
    //         lmin = evals[1];
    //         emax = evecs[0];
    //         emin = evecs[1];
    //     }
    // }

    field_t m_field;
    std::string m_name;
    pos_t m_default_h; 
};

template< typename Scalar_ >
struct planar_topology {
    typedef SmoothScalarField< Scalar_, 2 > field_t;
    typedef PWLScalarField< Scalar_, 2 > c0_field_t;
    typedef Scalar_ value_t;
    typedef typename field_t::pos_t pos_t;
    typedef typename field_t::gradient_t gradient_t;
    typedef typename field_t::hessian_t hessian_t;
    typedef typename field_t::bounds_t bounds_t;
    typedef typename field_t::solver_t solver_t;
    
    struct critical_point {
        pos_t x;
        int type; // 0: none, 1: saddle, 2: source, 3: sink, 
                  // 4: singular source, 5: singular sink, 6: 2nd order
        gradient_t emax, emin;
        value_t lmax, lmin;
        
        std::string as_string() const {
            switch (type) {
                case 0: return "invalid";
                case 1: return "saddle";
                case 2: return "source";
                case 3: return "sink";
                case 4: return "singular source";
                case 5: return "singular sink";
                case 6: return "2nd order critical point";
                default: return "error: unknown critical point type";
            }
        }
    };
    
    mutable std::list< std::pair< pos_t, gradient_t > > sampling_issues;
    mutable std::list< std::pair< pos_t, gradient_t > > ok_samples;
    mutable std::list< std::pair< pos_t, gradient_t > > current_samples;
    mutable std::list< std::pair< pos_t, gradient_t > > newton_steps;
    mutable std::list< std::pair< pos_t, gradient_t > > sampled_gradient;
    mutable std::list< pos_t > first_guesses;
    mutable std::list< std::list< pos_t > > search_path;
    
    planar_topology(const field_t& field) : m_field(field) {}
    
    const field_t& m_field;
        
    double edge_rotation(const pos_t& x0, const pos_t& x1, 
                         double dtheta, unsigned int max_depth,
                         unsigned int depth=0) const 
    {
        // needed for internal consistency
        const double large_angle = 0.75 * M_PI;
        
        gradient_t v0, v1, dummy;
        v0 = m_field.gradient(x0);
        v1 = m_field.gradient(x1);
        current_samples.push_back(std::make_pair(x0, v0));
        current_samples.push_back(std::make_pair(x1, v1));
        
        // easy cases first
        double theta = spurt::signed_angle(v0, v1);
        if (std::abs(theta) < std::max(large_angle, dtheta)) {
            return theta;
        }
        else if (depth > max_depth) {
            // use local linear extrapolation to figure out the rotation direction
            hessian_t h = m_field.hessian(x0);
            gradient_t dv0 = h * (x1 - x0);
            double outer_loc = v0[0] * dv0[1] - v0[1] * dv0[0];
            // verify that local rotation direction matches linear approximation along x0-x1
            double outer_edge = v0[0] * v1[1] - v0[1] * v1[0];
            if (outer_loc * outer_edge < 0) {
                std::ostringstream oss;
                oss << "Degenerate point encountered between " << flatten(x0) << " and " << flatten(x1) << ".";
                oss << "\nv0=" << flatten(v0) << ", v1=" << flatten(v1) << ", dv0=" << flatten(dv0);
                sampling_issues.insert(sampling_issues.end(), current_samples.begin(), current_samples.end());
                throw std::runtime_error(oss.str());
            }
            return theta;
        } 
        else {
            // use linear model to pick next sample
            // solve for v(x).v0 = 0
            double v0sq = spurt::vector::normsq(v0);
            double u = v0sq / (v0sq - spurt::vector::dot(v0, v1));
            pos_t x = (1 - u) * x0 + u * x1;
            return edge_rotation(
                       x0, x, dtheta, max_depth, depth + 1) +
                   edge_rotation(
                       x, x1, dtheta, max_depth, depth + 1);
        }
    }
    
    int index(const std::vector<pos_t>& pos) {
        double theta = 0;
        for (int i=0; i<4; ++i) {
            current_samples.clear();
            theta += edge_rotation(pos[i], pos[(i+1)%4], 0.01, 4);
            ok_samples.insert(ok_samples.end(), current_samples.begin(), current_samples.end());
        }
        double nb_rotations = 0.5 * theta / M_PI; 
        return int(nb_rotations);
    }
    
    struct box_constraint {
        bool inside(const pos_t& p) const {
            const pos_t& _min = m_box.first;
            const pos_t& _max = m_box.second;
            for (int dim=0; dim<2; ++dim) {
                if (p[dim] < _min[dim] || p[dim] > _max[dim]) return false;
            }
            return true;
        }
        
        pos_t clip(const pos_t& from, const pos_t& to) const {
            const pos_t& _min = m_box.first;
            const pos_t& _max = m_box.second;
            pos_t clipped(to);
            for (int dim=0; dim<2; ++dim) {
                if (clipped[dim] < _min[dim]) {
                    clipped[dim] = 0.05*from[dim] + 0.95*_min[dim];
                }
                else if (clipped[dim] > _max[dim]) {
                    clipped[dim] = 0.05*from[dim] + 0.95*_max[dim];
                }
            }
            return clipped;
        }
        
        box_constraint(const bounds_t& box) : m_box(box) {}
    
        pos_t operator()(const pos_t& from, const pos_t& to) const {
            assert(inside(from));
            if (inside(to)) {
                return to;
            }
            // std::cerr << "from = " << from << ", to = " << to << " (" << spurt::vector::norm(from-to) << ")" << std::endl;
            // std::cerr << "box = " << flatten(m_box.first) << " -> " << flatten(m_box.second) << std::endl;
            return clip(from, to);
        }
    
        double size() const {
            return nvis::norm(m_box.size());
        }
    
        bounds_t m_box;
    };
    
    bool lnsearch(pos_t& x, gradient_t& f, const pos_t& dd, double maxlength)
    {
        double lambda = 1.0;
        const double alpha = 1e-4;
    
        pos_t x0=x, xsave = x, f0=f, fsave = f;
        pos_t d = ( spurt::vector::norm(dd) > maxlength )
            ? dd * maxlength / spurt::vector::norm(dd) 
            : dd;
        value_t v0 = spurt::vector::norm(f);
        value_t vsave = v0;
        
        for (unsigned int i = 0; i < 7; ++i) {
            x = x0 + lambda * d;
            f = m_field.gradient(x);
            
            if (spurt::vector::norm(f) < (1 - alpha*lambda)*spurt::vector::norm(fsave)) {
                return true;
            }
        
            lambda *= 0.5;
        }
    
        return false;
    }

    bool newton_in_box(pos_t& x, const bounds_t& box,
                       double eps, size_t maxiter, bool verbose = false)
    {
        gradient_t d, f; // d is destination, f is rhs
        hessian_t H;
        pos_t delta = box.second - box.first;
        double maxlength = 0.5*delta.size();
        
        first_guesses.push_back(x);
        search_path.push_back( std::list< pos_t >() );
        std::list< pos_t >& polyline = search_path.back();
        
        box_constraint constraint(box);
    
        std::ostringstream os;
    
        pos_t best;
        double minnorm = std::numeric_limits<double>::max();
        bool must_improve = false;
        int nb_failed = 0;
    
        unsigned int k;
        try {
            f = m_field.gradient(x);
            double dinit = spurt::vector::norm(f);
            std::cout << "newton: seeding, norm(f(" << flatten(x) << ")) = " << dinit << '\n';
            minnorm = dinit;
            for (k = 0; k < maxiter; k++) {
                polyline.push_back(x);
                if (verbose) {
                    os << "newton: k = " << k << ", norm(f(" << flatten(x) << ")) = " << spurt::vector::norm(f)
                       << std::endl;
                    std::cerr << os.str();
                    os.clear();
                    os.str("");
                }
                double _norm = spurt::vector::norm(f);
                if (_norm < eps) return true;
                if (_norm < minnorm) {
                    minnorm = _norm;
                    best = x;
                    if (must_improve) nb_failed = 0;
                } else if (must_improve) {
                    ++nb_failed;
                    if (nb_failed > 5) {
                        break;
                    }
                }
            
                // determine local search direction
                H = m_field.hessian(x);
                d = solver_t::solve(H, H * x - f);
                
                // do a relaxation linesearch
                // (updates x and f)
                pos_t save(x);
                lnsearch(x, f, d - x, maxlength);
                // x = constraint(save, x);
                newton_steps.push_back(std::make_pair(save, save-x));
                sampled_gradient.push_back(std::make_pair(save, f));
            }
        
            if (k == maxiter) {
                if (verbose) {
                    os << "\t\t initial distance = " << dinit
                       << ", final distance = " << minnorm << ". failed.\n";
                    std::cerr << os.str();
                    os.clear();
                    os.str("");
                }
            }
        } catch (std::runtime_error& e) {
            if (verbose) {
                os << e.what() << std::endl;
                os << "\texception caught in Newton. current position is " << x << std::endl;
                std::cerr << os.str();
            }
            return false;
        }
    
        if (k == maxiter) {
            x = best;
        }
        return (minnorm < eps);
    }
    
    
    int resolve_critical_point(critical_point& cp, const std::vector<pos_t>& pos,
                               const pos_t& seed) {
        pos_t x = seed;
        bounds_t box(pos[0], pos[2]);
        if (newton_in_box(x, box, 1.0e-6, 50, true)) {
            cp.x = x;
            m_field.hessian_eigen(x, cp.lmax, cp.lmin, cp.emax, cp.emin);
            if (cp.lmin * cp.lmax < 0) cp.type = 1;
            else if (cp.lmin > 0) cp.type = 2;
            else if (cp.lmax < 0) cp.type = 3;
            else if (cp.lmin == 0 && cp.lmax > 0) cp.type = 4;
            else if (cp.lmax == 0 && cp.lmin < 0) cp.type = 5;
            else cp.type = 6;
            return cp.type;
        }
        return 0;
    }
    
    void find_critical_points(std::vector<critical_point>& cpts, 
                              size_t nx=10, size_t ny=10) {
        cpts.clear();
        bounds_t bounds = m_field.bounds();
        pos_t spacing = bounds.second - bounds.first;
        const pos_t& origin = bounds.first;
        spacing[0] /= nx;
        spacing[1] /= ny;
        std::vector< pos_t > pos(4);
        int nfound = 0;
        for (int j=0; j<ny; ++j) {
            double ymin = origin[1] + j * spacing[1];
            double ymax = origin[1] + (j+1) * spacing[1];
            for (int i=0; i<nx; ++i) {
                double xmin = origin[0] + i * spacing[0];
                double xmax = origin[0] + (i + 1) * spacing[0];
                pos[0] = pos_t(xmin, ymin);
                pos[1] = pos_t(xmax, ymin);
                pos[2] = pos_t(xmax, ymax);
                pos[3] = pos_t(xmin, ymax);
                bool failed = false;
                try {
                    int indx = index(pos); 
                    if (indx != 0) {
                        std::cout << "critical point with index " << indx << " found in cell (" << i << ", " << j << ")\n";
                        ++nfound;
                        critical_point cp;
                        pos_t seed = 0.5*(pos[0] + pos[2]);
                        if (!resolve_critical_point(cp, pos, seed)) {
                            std::cerr << "Unable to find critical point\n";
                            failed = true;
                        }
                        
                        int res = 4;
                        for (int n=0; n<3 && failed; ++n, res*=2) {
                            std::cerr << "Trying again with better seed (res=" << res << ")\n";
                            pos_t dx = 1./res*(pos[2]-pos[0]);
                            pos_t best = pos[0];
                            value_t bestv = std::numeric_limits<value_t>::max();
                            for (int l=0; l<(res+1); ++l) {
                                for (int k=0; k<(res+1); ++k) {
                                    pos_t p = pos[0] + pos_t(l*dx[0], k*dx[1]);
                                    value_t v = spurt::vector::norm(m_field.gradient(p));
                                    if (v < bestv) {
                                        bestv = v;
                                        best = p;
                                        std::cerr << "new min norm found at " << flatten(p) << ", norm=" << v << '\n';
                                    }
                                }
                            }
                            if (resolve_critical_point(cp, pos, best)) {
                                failed = false;
                            }
                            else {
                                std::cerr << "unable to find critical point despite repeated attempts\n";
                            }
                        }
                        
                        if (!failed) {
                            std::cout << "critical point of type " << cp.as_string() << " found!\n";
                            std::cout << "\tposition: " << flatten(cp.x) << '\n';
                            std::cout << "\tlmax: " << cp.lmax << '\n';
                            std::cout << "\tlmin: " << cp.lmin << '\n';
                            std::cout << "\temax: " << flatten(cp.emax) << '\n';
                            std::cout << "\temin: " << flatten(cp.emin) << '\n';
                            gradient_t f = m_field.gradient(cp.x);
                            std::cout << "\tvalue: " << flatten(f) << '\n';
                            std::cout << "\tnorm: " << spurt::vector::norm(f) << '\n';
                            cpts.push_back(cp);
                        }
                    }
                }
                catch(std::exception& e) {
                    // std::cerr << "exception caught while computing index in cell (" << i << ", " << j << ")\n";
                    // std::cerr << e.what() << std::endl;
                    // std::cerr << "skipping cell\n";
                    // std::cerr << sampling_issues.size() << " samples currently in list\n";
                }
            }
        }
        
        std::cout << "\n\n" << nfound << " critical points found\n";
        if (!sampling_issues.empty()) {
            std::list< pos_t > pos;
            std::list< gradient_t > vecs;
            std::for_each(sampling_issues.begin(), sampling_issues.end(), 
            [&](const std::pair<pos_t, gradient_t>& p)
                {
                   pos.push_back(p.first);
                   vecs.push_back(p.second); 
                });
            VTK_SMART(vtkPolyData) poly = vtk_utils::make_points(pos);
            poly = vtk_utils::add_vectors(poly, vecs);
            vtk_utils::saveVTK(poly, "sampling_issues.vtp");
        }
        
        if (!ok_samples.empty()) {
            std::list< pos_t > pos;
            std::list< gradient_t > vecs;
            std::for_each(ok_samples.begin(), ok_samples.end(), 
            [&](const std::pair<pos_t, gradient_t>& p)
                {
                   pos.push_back(p.first);
                   vecs.push_back(p.second); 
                });
            VTK_SMART(vtkPolyData) poly = vtk_utils::make_points(pos);
            poly = vtk_utils::add_vectors(poly, vecs);
            vtk_utils::saveVTK(poly, "other_samples.vtp");
        }
        
        if (!newton_steps.empty()) {
            std::list< pos_t > pos;
            std::list< gradient_t > vecs;
            std::for_each(newton_steps.begin(), newton_steps.end(), 
            [&](const std::pair<pos_t, gradient_t>& p)
                {
                   pos.push_back(p.first);
                   vecs.push_back(p.second); 
                });
            VTK_SMART(vtkPolyData) poly = vtk_utils::make_points(pos);
            poly = vtk_utils::add_vectors(poly, vecs);
            vtk_utils::saveVTK(poly, "newtonsteps.vtp");
        }
        
        if (!sampled_gradient.empty()) {
            std::list< pos_t > pos;
            std::list< gradient_t > vecs;
            std::for_each(sampled_gradient.begin(), sampled_gradient.end(), 
            [&](const std::pair<pos_t, gradient_t>& p)
                {
                   pos.push_back(p.first);
                   vecs.push_back(p.second); 
                });
            VTK_SMART(vtkPolyData) poly = vtk_utils::make_points(pos);
            poly = vtk_utils::add_vectors(poly, vecs);
            vtk_utils::saveVTK(poly, "gradients.vtp");
        }
        
        if (!first_guesses.empty()) {
            VTK_SMART(vtkPolyData) poly = vtk_utils::make_points(first_guesses);
            poly = vtk_utils::add_vertices(poly);
            vtk_utils::saveVTK(poly, "seeds.vtp");
        }

        if (!search_path.empty()) {
            VTK_SMART(vtkPolyData) poly = vtk_utils::make_polylines(search_path);
            vtk_utils::saveVTK(poly, "search.vtp");
        }
                   
    }
};


} // scalar 
} // topology
} // spurt





#endif