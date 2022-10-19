// #include "extractor.hpp"
#include "math/matrix.hpp"
#include <sstream>
#include <string>
#include "pvo.hpp"
#include "measure_wrapper.hpp"
#include "face.hpp"

#define QUICK_TEST_3

using namespace xavier;
using namespace xavier::crease;

const double DET_EPSILON = 1.0e-6;
const double inside_eps = 0.01;

bool __verbose = false;
nvis::vec3 actual_pos;
unsigned int __max_depth;

double current_face_avg_norm;

unsigned int xavier::crease::nb_pvo = 0;

inline double dist(double v)
{
    return std::max(std::max(-v, 0.), std::max(v - 1, 0.));
}

// ---------------------------------------------------------------------------

// Stetten's average direction
nvis::vec3 xavier::crease::average(nvis::vec3 dirs[4])
{
    double mat[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    // mat = sum_i( dirs[i] * dirs[i]^T )
    for (unsigned int l = 0 ; l < 4 ; l++) {
        for (unsigned int i = 0 ; i < 3 ; i++)
            for (unsigned int j = i ; j < 3 ; j++) {
                mat[i+3*j] += 0.25 * dirs[l][i] * dirs[l][j];
            }
    }
    mat[3] = mat[1];
    mat[6] = mat[2];
    mat[7] = mat[5];
    
    double eval[3];
    double evec[9];
    int nbroots = ell_3m_eigensolve_d(eval, evec, mat, 1);
    
    nvis::vec3 avg(evec[0], evec[1], evec[2]);
    
    return avg;
}

// ---------------------------------------------------------------------------

face_type* current_face;

std::ostream& operator<<(std::ostream& out, const std::vector< int >& refs)
{
    for (unsigned int i = 0 ; i < refs.size() ; ++i) {
        out << refs[i];
        if (i < refs.size() - 1) {
            out << "-";
        }
    }
    return out;
}

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

void subdivide_face(std::vector< face_type >& out, const face_type& face, bool interpolate)
{
    out.clear();
    nvis::vec3 p[9], g[9], Hg[9];
    for (unsigned int i = 0 ; i < 9 ; ++i) {
        p[i] = position(i, face);
        g[i] = gradient(i, face, interpolate);
        Hg[i] = Hgradient(i, face, interpolate);
    }
    
    for (unsigned int i = 0 ; i < 4 ; i++) { // loop over subquads
        out.push_back(face_type());
        face_type& f = out.back();
        std::copy(face.reference.begin(), face.reference.end(), std::back_inserter(f.reference));
        f.reference.push_back(i);
        f.depth = face.depth + 1;
        
        for (unsigned int j = 0 ; j < 4 ; j++) { // loop over vertices
            f.p[j] = p[xavier::crease::indices[i][j]];
            f.g[j] = g[xavier::crease::indices[i][j]];
            f.Hg[j] = Hg[xavier::crease::indices[i][j]];
        }
        f.set_basis();
    }
}

// ---------------------------------------------------------------------------

void xavier::crease::refine_face(face_type& out, const face_type& face,
                                 const nvis::vec3& q, double h)
{
    // compute local coordinates
    nvis::vec2 x = local_coord(face, q);
    double xmin = std::max(0., x[0] - h);
    double xmax = std::min(1., x[0] + h);
    double ymin = std::max(0., x[1] - h);
    double ymax = std::min(1., x[1] + h);
    out = face_type(face(xmin, ymin), face(xmax, ymin),
                    face(xmax, ymax), face(xmin, ymax));
    std::copy(face.reference.begin(), face.reference.end(), std::back_inserter(out.reference));
    out.reference.push_back(8);
    
    out.depth = face.depth;
    // we do not increment the depth when we refine around a PVO solution
}

// ---------------------------------------------------------------------------

bool error3(const nvis::vec3& ref, const nvis::vec3& val, bool normalize = true)
{
    using namespace xavier::crease;
    const double min_dot = cos(max_int_error * M_PI / 2.);
    
    double mr = 1, mv = 1;
    if (normalize) {
        mr = nvis::norm(ref);
        mv = nvis::norm(val);
    }
    return (mr > gradient_eps &&
            nvis::inner(ref, val) / (mr*mv) < min_dot);
}

bool error2(const face_type& f, double ref_g_norm, double ref_Hg_norm)
{
    using namespace xavier::crease;
    
    double gm = std::max(gradient_eps, max_int_error * ref_g_norm);
    double Hgm = std::max(gradient_eps, max_int_error * ref_Hg_norm);
    nvis::vec3 g[2], Hg[2];
    g[0] = gradient(8, f, false);
    g[1] = gradient(8, f, true);
    Hg[0] = Hgradient(8, f, false);
    Hg[1] = Hgradient(8, f, true);
    
    double dg = nvis::norm(g[0] - g[1]);
    double gradmag = nvis::norm(g[0]);
    if (gradmag > gradient_eps && // failed vanishing gradient test
            dg > gm && // failed distance test
            nvis::inner(g[0], g[1]) < 1. - max_int_error) { // failed alignment test
        if (__verbose) {
            std::cout << "In face " << f.reference
                      << ": gradient error at midpoint = " << dg << " exceeds threshold (" << gm << ")" << std::endl
                      << "mean gradient = " << g[0] << " (" << nvis::norm(g[0]) << "), values are: " << std::flush;
            for (unsigned int i = 0 ; i < 4 ; ++i) {
                std::cout << i << ": " << f.g[i] << " (" << nvis::norm(f.g[i]) << "), ";
            }
            std::cout << std::endl;
        }
        return true;
    }
    
    double dHg = nvis::norm(Hg[0] - Hg[1]);
    if (dHg > Hgm && // failed distance test
            nvis::inner(Hg[0], Hg[1]) < 1. - max_int_error) { // failed alignment test
        if (__verbose) {
            std::cout << "In face " << f.reference
                      << ": Hgradient error at midpoint = " << dHg << " exceeds threshold (" << Hgm << ")" << std::endl
                      << "mean Hgradient = " << Hg[0] << " (" << nvis::norm(Hg[0]) << "), values are: " << std::flush;
            for (unsigned int i = 0 ; i < 4 ; ++i) {
                std::cout << i << ": " << f.Hg[i] << " (" << nvis::norm(f.Hg[i]) << "), ";
            }
            std::cout << std::endl;
        }
        return true;
    }
    
    return false;
}

// ---------------------------------------------------------------------------

bool xavier::crease::par_vec_op(std::vector< nvis::vec3 >& beta, const xavier::mat3& V,
                                const xavier::mat3& W)
{
    beta.clear();
    
    using namespace xavier;
    using namespace crease;
    
    // check if W can be inverted
    xavier::mat3 M;
    
    
    double detW = fabs(det(W));
    double detV = fabs(det(V));
    if (detW >= detV) {
        if (detW == 0) {
            return false;
        }
        M = invert(W) * V;
    } else if (detV > detW) {
        M = invert(V) * W;
    } else {
        return false;
    }
    
    // compute eigenvectors of resulting matrix
    std::vector< double > evals;
    std::vector< nvis::vec3 > evecs;
    M.eigensystem(evals, evecs); // returns only real eigenvalues
    
    // checking the results
    if (!evals.size()) {
        if (__verbose) {
            std::cout << "no real eigenvalue to PVO eigensystem" << std::endl;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < evals.size() ; i++) {
        if (std::isnan(evals[i])) {
            if (__verbose) {
                std::cout << "NaN eigenvalue to PVO eigensystem" << std::endl;
            }
            return false;
        }
    }
    
    // looking for valid solutions
    bool found = false;
    for (unsigned int i = 0 ; i < evecs.size() ; i++) {
        nvis::vec3 _ev = evecs[i];
        double scal = _ev[2];
        if (scal == 0) {
            continue;
        }
        
        _ev *= 1 / scal;
        if (_ev[0] > -inside_eps && _ev[1] > -inside_eps && _ev[0] + _ev[1] < 1 + inside_eps) {
            found = true;
            beta.push_back(nvis::vec3(1 - _ev[0] - _ev[1], _ev[0], _ev[1]));
            
            // NB: we do not verify that the found location does indeed
            // correspond to the expected kind of crease point.
            // This must be done afterwards
            
            if (__verbose) {
                std::cout << "PVO found valid solution #" << beta.size() << ": " << beta.back()
                          << std::endl;
            }
        } else if (fixing_voxel && __verbose) {
            std::cout << "wrong coordinates in PVO: "
                      << 1 - _ev[0] - _ev[1] << ", " << _ev[0] << ", " << _ev[1] << std::endl;
        }
    }
    
    return found;
}

// ---------------------------------------------------------------------------

bool xavier::crease::linear_parallel_operator(std::vector< nvis::vec3 >& b,
        const nvis::vec3 g[3],
        const nvis::vec3 ev[3])
{
    // parallel vector operator method for linear faces
    // use barycentric coordinates: #0->(0,0), #1->(1,0), #2->(0,1)
    // solve for ev_i = V*(s,t,1)^T, g_i = W*(s,t,1)^T
    xavier::mat3 V, W;
    for (unsigned int i = 0 ; i < 3 ; i++) {
        V(i, 2) = ev[0][i];
        V(i, 0) = ev[1][i] - ev[0][i];
        V(i, 1) = ev[2][i] - ev[0][i];
        
        W(i, 2) = g[0][i];
        W(i, 0) = g[1][i] - g[0][i];
        W(i, 1) = g[2][i] - g[0][i];
    }
    
    if (__verbose) {
        std::cout << "linear_parallel_vector_operator:" << std::endl
                  << "input gradients: g[0] = " << g[0] << ", g[1] = " << g[1] << ", g[2] = " << g[2] << std::endl
                  << "input eigenvectors: e[0] = " << ev[0] << ", e[1] = " << ev[1] << ", e[2] = " << ev[2] << std::endl;
    }
    
    return par_vec_op(b, V, W);
}

// ---------------------------------------------------------------------------

bool check_face_for_crossings(const face_type& face)
{
    if (fixing_voxel) {
        return true;    // no smarty pants testing in second pass
    }
    
#ifdef QUICK_TEST_1
    nvis::vec2 g_tan[4], Hg_tan[4];
    
    unsigned int nb_crossings = 0;
    
    for (unsigned int d = 0 ; d < 3 ; d++) {
        unsigned int d1, d2;
        d1 = d;
        d2 = (d + 1) % 3;
        
        for (unsigned int i = 0 ; i < 4 ; i++) {
            g_tan[i][0] = face.g[i][d1];
            g_tan[i][1] = face.g[i][d2];
            Hg_tan[i][1] = face.Hg[i][d2];
            Hg_tan[i][0] = face.Hg[i][d1];
        }
        
        double cross[4];
        for (unsigned int i = 0 ; i < 4 ; i++) {
            cross[i] = g_tan[i][0] * Hg_tan[i][1] - g_tan[i][1] * Hg_tan[i][0];
        }
        
        unsigned int nb_zeros = 0;
        
        for (unsigned int i = 0 ; i < 4 ; i++) {
            unsigned int j = (i + 1) % 4;
            if (cross[i]*cross[j] < 0) {
                ++nb_zeros;
            }
        }
        
        if (nb_zeros == 2) {
            ++nb_crossings;
        }
    }
    
    return (nb_crossings == 3);
#endif
    
#ifdef QUICK_TEST_2
    
    // first compute average gradient and eigenvector direction
    std::vector< nvis::vec3 > basis(3);
    crease::PCA(basis, xprod);
    
    // express cross products in the basis formed by the 2 main eigendirections
    nvis::vec2 x2d[4];
    for (unsigned int i = 0 ; i < 4 ; i++) {
        for (unsigned int j = 0 ; j < 2 ; j++) {
            x2d[i][j] = nvis::inner(xprod[i], basis[j]);
        }
    }
    
    // compute the number of zero crossings
    unsigned int nbcross[2] = {0, 0};
    for (unsigned int i = 0; i < 4 ; i++) {
        unsigned int ii = (i + 1) % 4;
        for (unsigned int j = 0 ; j < 2 ; j++) {
            if (x2d[i][j]*x2d[ii][j] <= 0) {
                nbcross[j]++;
            }
        }
    }
    
    // we need at least a zero crossing for each dimension
    return (nbcross[0]*nbcross[1] > 0);
    
#endif
    
#ifdef QUICK_TEST_3
    
    if (true || !fixing_voxel) {
        nvis::vec3 xprod[4];
        for (unsigned int i = 0 ; i < 4 ; i++) {
            xprod[i] = nvis::cross(face.g[i], face.Hg[i]);
        }
        
        int nbcross[3] = {0, 0, 0};
        for (unsigned int i = 0 ; i < 4 ; ++i) {
            unsigned int ii = (i + 1) % 4;
            for (unsigned int j = 0 ; j < 3 ; ++j) {
                if (xprod[i][j]*xprod[ii][j] <= 0) {
                    nbcross[j]++;
                }
            }
        }
        
        if (true || xavier::crease::speedup) {
            return (nbcross[0]*nbcross[1] + nbcross[0]*nbcross[2]  + nbcross[1]*nbcross[2] > 0);
        } else {
            return (nbcross[0] + nbcross[1] + nbcross[2] > 0);
        }
    } else {
        std::vector< nvis::vec3 > xprod(4);
        for (unsigned int i = 0 ; i < 4 ; i++) {
            xprod[i] = nvis::cross(face.g[i], face.Hg[i]);
        }
        
        nvis::vec3 main_dir = crease::average(xprod);
        double prod[4];
        for (unsigned int i = 0 ; i < 4 ; i++) {
            prod[i] = nvis::inner(xprod[i], main_dir);
        }
        for (unsigned int i = 0 ; i < 4 ; ++i) {
            unsigned int j = (i + 1) % 4;
            if (prod[i]*prod[j] <= 0) {
                return true;
            }
        }
        
        return false;
    }
    
#endif
    
#ifdef QUICK_TEST_4
    
    std::vector< nvis::vec3 > xprod(4);
    for (unsigned int i = 0 ; i < 4 ; i++) {
        xprod[i] = nvis::cross(face.g[i], face.Hg[i]);
    }
    
    nvis::vec3 main_dir = crease::average(xprod);
    double prod[4];
    for (unsigned int i = 0 ; i < 4 ; i++) {
        prod[i] = nvis::inner(xprod[i], main_dir);
    }
    for (unsigned int i = 0 ; i < 4 ; ++i) {
        unsigned int j = (i + 1) % 4;
        if (prod[i]*prod[j] <= 0) {
            return true;
        }
    }
    
#endif
    
}

// ---------------------------------------------------------------------------

//  1: sucessful
//  0: nothing found
bool process_face_PVO(std::vector< nvis::vec3 >& xing,
                      const face_type& face, const unsigned int depth)
{
    using namespace xavier::crease;
    
    ++nb_pvo;
    
    // std::cout << "process_face: d=" << depth << std::endl;
    // std::cout << "face pos: " << p0 << " " << p1 << " " << p2 << " " << p3 << std::endl;
    bool found = false;
    
    // bool verbose = !xavier::crease::apply_filter;
    bool verbose = __verbose;
    
    if (verbose) {
        std::cout << "PVO in face " << face.reference << std::endl
                  << "processing face: " << face.p[0] << ", " << face.p[1]
                  << ", " << face.p[2] << ", " << face.p[3] << std::endl;
    }
    
    // add current quadrilateral to the list
    if (xavier::crease::display_debug_info) {
        crease::current_vertices.push_back(face.p[0]);
        crease::current_vertices.push_back(face.p[1]);
        crease::current_vertices.push_back(face.p[2]);
        crease::current_vertices.push_back(face.p[3]);
    }
    // good approximation quality achieved through bilinear interpolation
    // extract parallel location using PVO
    unsigned int triangles[2][3] = { { 0, 1, 2 }, { 0, 2, 3 } };
    nvis::vec3 g[3], Hg[3];
    
    for (unsigned int t = 0 ; t < 2 ; t++) {
        std::vector< nvis::vec3 > bs;
        
        g[0] = face.g[triangles[t][0]];
        g[1] = face.g[triangles[t][1]];
        g[2] = face.g[triangles[t][2]];
        
        Hg[0] = face.Hg[triangles[t][0]];
        Hg[1] = face.Hg[triangles[t][1]];
        Hg[2] = face.Hg[triangles[t][2]];
        
        if (linear_parallel_operator(bs, g, Hg)) {
        
            // if (xavier::crease::display_debug_info) {
            //  xavier::crease::pvo_faces.push_back(face.p[0]);
            //  xavier::crease::pvo_faces.push_back(face.p[1]);
            //  xavier::crease::pvo_faces.push_back(face.p[2]);
            //  xavier::crease::pvo_faces.push_back(face.p[3]);
            // }
            
            for (unsigned int n = 0 ; n < bs.size() ; ++n) {
                nvis::vec3& beta = bs[n];
                if (__verbose)
                    std::cout << "d: " << depth
                              << "... PVO successful in triangle #" << t << "... " << std::endl;
                              
                nvis::vec3 q = beta[0] * face.p[triangles[t][0]] +
                               beta[1] * face.p[triangles[t][1]] +
                               beta[2] * face.p[triangles[t][2]];
                xing.push_back(q);
                found = true;
                
                // debugging
                if (__verbose) {
                    std::cout << "*************************************" << std::endl;
                    std::cout << "PVO found a solution at " << q << " in triangle "
                              << triangles[t][0] << ", " << triangles[t][1] << ", " << triangles[t][2] << std::endl;
                    nvis::vec3 _g, _e;
                    _g = crease::the_wrapper->gradient(q);
                    _e = crease::the_wrapper->eigenvector(q, 0);
                    std::cout << "ground truth:" << std::endl;
                    std::cout << "value is " << std::setprecision(20) << crease::the_wrapper->value(q) << std::endl;
                    std::cout << "gradient is " << std::setprecision(8) << _g << std::endl;
                    std::cout << "eigenvector is " << _e << std::endl;
                    std::cout << "cross product norm is " << nvis::norm(nvis::cross(_g, _e)) << std::endl;
                    std::cout << "local linear approximation:" << std::endl;
                    std::cout << "barycentric coordinates are " << beta << std::endl;
                    nvis::vec3 grad = beta[0] * g[0] + beta[1] * g[1] + beta[2] * g[2];
                    nvis::vec3 Hgrad = beta[0] * Hg[0] + beta[1] * Hg[1] + beta[2] * Hg[2];
                    double gm = nvis::norm(grad);
                    double Hgm = nvis::norm(Hgrad);
                    std::cout << "gradient = " << grad << " (norm=" << gm << ")" << std::endl
                              << "the three gradient values are: #0: " << g[0] << ", #1: " << g[1] << ", #2:" << g[2] << std::endl
                              << "the three Hgradient values are: #0: " << Hg[0] << ", #1: " << Hg[1] << ", #2:" << Hg[2] << std::endl
                              << "dot product = " << fabs(nvis::inner(grad, Hgrad)) / gm << std::endl
                              << "cross product = " << nvis::norm(nvis::cross(Hgrad, grad)) << std::endl;
                    std::cout << "*************************************" << std::endl;
                }
            }
        } else if (__verbose) {
            std::cout << "PVO was unsuccessful in triangle #" << t << " of face " << face.reference << std::endl;
        }
    }
    
    return found;
}

unsigned int bad_counter = 0;
int xavier::crease::
search_face_PVO::operator()(std::vector< nvis::vec3 >& xing,
                            const nvis::vec3& p0, const nvis::vec3& p1,
                            const nvis::vec3& p2, const nvis::vec3& p3,
                            unsigned int maxdepth, bool& something_found) const
{
    something_found = false;
    
    // std::cerr << "search_face_PVO: " << p0 << " - " << p2 << "\n";
    
    __max_depth = maxdepth;
    if (fixing_voxel && display_debug_info) {
        __verbose = true;
    } else {
        __verbose = false;
    }
    
    bool found = false;
    unsigned int depth = 0;
    xing.clear();
    current_vertices.clear();
    std::vector< nvis::vec3 > cur_xing;
    std::vector< face_type > cur, next, sub_approx, sub_measure;
    face_type ref_face(p0, p1, p2, p3);
    cur.push_back(ref_face);
    bool say = false;
    
    double ref_grad_mag = ref_face.average_norm(ref_face.g);
    double ref_Hgrad_mag = ref_face.average_norm(ref_face.Hg);
    
    unsigned int nbsq = 0;
    unsigned int bad_counter_aux = 0;
    
    // we allow only one crease point per face since we will not be able
    // to make sense of the mess if we find multiple solutions
    while (cur.size() && !found) {
    
        for (unsigned i = 0 ; i < cur.size() && !found ; i++) {
        
            const face_type& fp = cur[i];
            current_face = &cur[i];
            
            unsigned int depth = fp.depth;
            if (depth > maxdepth) {
                continue;    // should not happen!
            }
            
            bool bad = false;
            if (depth < maxdepth && // can still subdivide
                    error2(fp, ref_grad_mag, ref_Hgrad_mag)) {
//              (error3(gref, gapprox) || error3(eref, eapprox))) {
                subdivide_face(sub_measure, fp, false);
                bad = true;
            }
            
            if (bad) {
                // poor approximation quality - subdivision required
                if (__verbose) {
                    std::cout << "poor approximation quality in ";
                    for (unsigned int l = 0 ; l < fp.reference.size() ; ++l) {
                        std::cout << fp.reference[l] << "-";
                    }
                    std::cout << ": subdividing at depth " << depth << std::endl;
                }
                std::copy(sub_measure.begin(), sub_measure.end(),
                          std::back_inserter(next));
            } else {
                if (!check_face_for_crossings(fp)) {
                    continue;
                }
                ++nbsq;
                cur_xing.clear();
                if (process_face_PVO(cur_xing, fp, depth)) {
                    if (__verbose) {
                        std::cout << "looping over found PVO solutions (" << cur_xing.size() << ")" << std::endl;
                        if (cur_xing.size() > 1) {
                            std::cout << "we obtained several solutions!" << std::endl;
                        }
                    }
                    for (unsigned int n = 0 ; n < cur_xing.size() ; n++) {
                    
                        if (__verbose) {
                            std::cout << "PVO solution #" << n << std::endl;
                        }
                        
                        nvis::vec3 q = cur_xing[n];
                        // measures
                        double val = xavier::crease::the_wrapper->value(q);
                        nvis::vec3 g = xavier::crease::the_wrapper->gradient(q);
                        nvis::vec3 Hg = xavier::crease::the_wrapper->Hgradient(q);
                        nvis::vec3 e = xavier::crease::the_wrapper->eigenvector(q, xavier::crease::is_ridge ? 0 : 2);
                        double str = xavier::crease::the_wrapper->eigenvalue(q, 1);
                        // quality metrics
                        double gm = nvis::norm(g);
                        double Hgm = nvis::norm(g);
                        double dot = fabs(nvis::inner(g, Hg)) / gm / Hgm;
                        double cross = nvis::norm(nvis::cross(g, Hg));
                        double avggm = fp.average_norm(fp.g);
                        double avgHgm = fp.average_norm(fp.Hg);
                        
                        // special case for FA ridge lines
                        if (crease_kind == 1 && val > 1.01) {
                            something_found = true;
                            return false;
                        }
                        
                        bool tmp_found = false;
                        while (true) {
                            // did we find an extremum?
                            if (gm < std::max(xavier::crease::gradient_eps_rel*avggm,
                                              xavier::crease::gradient_eps)) {
                                tmp_found = true;
                                break;
                            }
                            
                            // check that we have a valid solution of the PVO
                            if (dot > (1 - xavier::crease::max_align_error)) {
                                // is this the eigenvector we were looking for
                                if (fabs(nvis::inner(g, e)) > 1 - xavier::crease::max_align_error) {
                                    tmp_found = true;
                                    break;
                                } else {
                                    tmp_found = false;
                                    break;
                                }
                            } else {
                                // linear approximation of poor quality
                                if (__verbose) {
                                    std::cout << "inaccurate PVO solution - refine around found position" << std::endl;
                                }
                                face_type tmp;
                                refine_face(tmp, fp, q);
                                next.push_back(tmp);
                                if (__verbose) {
                                    std::cout << "\t inaccurate PVO solution at " << q << " (" << dot << ", "
                                              << gm << ", " << 100*gm / avggm << "%): refining around it at depth "
                                              << depth << std::endl;
                                    say = true;
                                }
                                tmp_found = false;
                                break;
                            }
                            break;
                        }
                        
                        if (!tmp_found) {
                            continue;
                        }
                        
                        // 2nd: check crease strength
                        if (!is_ok(val, str)) {
                            // the solution is valid but irrelevant
                            if (__verbose)
                                std::cout << "\t PVO valid solution at " << q << " (dot=" << dot
                                          << ", gm=" << gm
                                          << ", val=" << val << ") is irrelevant: value is " << val << "("
                                          << crease::value_threshold_select << ") and strength is " << str
                                          << "(" << crease::strength_threshold_select << "), depth is " << depth << std::endl;
                                          
                            // make sure that we did not get ahead of ourselves with our shaky heuristics
                            if (gm > xavier::crease::gradient_eps &&
                                    dot < (1 - xavier::crease::max_align_error)) {
                                // this is still a pretty bad approximation - let us try to do better
                                // before we give up
                                face_type tmp;
                                refine_face(tmp, fp, q);
                                next.push_back(tmp);
                                if (__verbose) {
                                    std::cout << "\t inaccurate PVO solution at " << q << " (" << dot << ", "
                                              << gm << ", " << 100*gm / avggm << "%): refining around it at depth "
                                              << depth << std::endl;
                                    say = true;
                                }
                            } else {
                                something_found = true;
                            }
                        }
                        // 3rd: check that the position we found lies in the face
                        else if (!inside(ref_face, q)) {
                            // the solution is valid, relevant, but outside of the cell
                            if (__verbose) {
                                std::cout << "\t PVO solution is good but outside of face" << std::endl;
                            }
                        }
                        // this is the real thing
                        else {
                            std::cout << "found crease point around at " << q << " and depth " << depth
                                      << std::endl << "value = " << xavier::crease::the_wrapper->value(q)
                                      << std::endl << "<g,e>/|g| = " << fabs(nvis::inner(g, e)) << std::endl
                                      << "gm = " << nvis::norm(g) << std::endl
                                      << "strength = " << str << std::endl;
                                      
                            something_found = true;
                            
                            xing.push_back(q);
                            found = true;
                            /*
                            // an aggressive strategy is in place that interrupts the search after
                            // the first point found
                            // return 1;
                            */
                            break;
                        }
                    }
                }
            }
        }
        std::swap(cur, next);
        next.clear();
    }
    
    if (found)
        std::cout << "final depth = " << depth
                  << " (max was: " << maxdepth << "), "
                  << nbsq << " quads processed" << std::endl;
                  
    bad_counter = bad_counter_aux;
    
    return (found ? 1 : 0);
}





















