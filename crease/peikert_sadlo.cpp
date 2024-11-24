// #include <crease/peikert_sadlo.hpp>
#include <vtk/vtk_utils.hpp>
#include <vtk/vtk_interpolator.hpp>
#include <boundary-aware-rectgrid/boundaryAwareRectGrid.h>

#include <misc/option_parse.hpp>
#include <misc/sort.hpp>
#include <Eigen/Eigen>

#include <misc/progress.hpp>
#include <math/types.hpp>

using namespace spurt;

typedef double value_type;
typedef vec3   position_type;
typedef vec3   vector_type;
typedef cvec3  complex_vector_type;
typedef cmat3  complex_matrix_type;
typedef mat3   matrix_type;
typedef vec2   vector2d_type;
typedef vec2   position2d_type;

constexpr value_type invalid_value = std::numeric_limits<value_type>::min();
inline bool is_bad_value(value_type v) {
    return v == invalid_value;
}

typedef lvec2 edge_index_type;
typedef lvec3 triangle_index_type;
typedef lvec4 quad_index_type;
typedef lvec4 face_index_type;

typedef vtk_utils::interpolator<vtkUnstructuredGrid, double, 3, vec3, mat3> 
    unstructured_interpolator_type;
typedef vtk_utils::interpolator<vtkStructuredGrid, double, 3, vec3, mat3> 
    curvilinear_interpolator_type;
typedef vtk_utils::interpolator<vtkRectilinearGrid, double, 3, vec3, mat3> 
    rectilinear_interpolator_type;


std::string name_in, name_out, name_surface;
double eps;
double strength_threshold;
double value_threshold;
bool verbose;
int niter;
bool use_newton;
bool do_ridges;
bool do_valleys;
bool do_surfaces;
bool do_lines;
bool do_ridge_lines, do_ridge_surfaces, do_valley_lines, do_valley_surfaces;
bool do_omega, do_l2, do_l, do_h;

void initialize(int argc, char* argv[]) {
    namespace scl = spurt::command_line;

    scl::option_traits
            required(true, false, "Required Options"),
            optional(false, false, "Optional Group");
    scl::option_parser parser(argv[0],
            "Extract ridges of scalar field using Peikert and Sadlo's method");

    verbose = false;
    eps = 1.0e-6;
    niter = 20;
    use_newton = true;
    do_ridges = true;
    do_valleys = true;
    do_surfaces = true;
    do_lines = true;
    value_threshold = 0;
    strength_threshold = 0;
    name_out = "creases.vtp";
    do_omega = true;
    do_l2 = true;
    do_l = true;
    do_h = true;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required);
        parser.add_value("output", name_out, name_out, "Output base name", optional);
        parser.add_value("ridges", do_ridges, do_ridges, "Extract ridge creases", optional);
        parser.add_value("valleys", do_valleys, do_valleys, "Extract valley creases", optional);
        parser.add_value("surfaces", do_surfaces, do_surfaces, "Extract surface creases", optional);
        parser.add_value("lines", do_lines, do_lines, "Extract line creases", optional);
        parser.add_value("vorticity", do_omega, do_omega, "Compute vorticity", optional);
        parser.add_value("lambda2", do_l2, do_l2, "Compute lambda2", optional);
        parser.add_value("lamb", do_l, do_l, "Compute lamb vector", optional);
        parser.add_value("helicity", do_h, do_h, "Compute helicity", optional);
        parser.add_value("surface", name_surface, name_surface, "Reampling surface", optional);
        parser.add_value("value", value_threshold, value_threshold, "Value threshold");
        parser.add_value("strength", strength_threshold, strength_threshold, "Crease strength threshold");
        parser.add_value("eps", eps, eps, "Integration precision", optional);
        parser.add_value("newton", use_newton, use_newton, "Use Newton search", optional);
        parser.add_value("niter", niter, niter, "Number of Newton iterations", optional);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);

        parser.parse(argc, const_cast<const char**>(argv));

        do_ridge_lines = do_lines && do_ridges;
        do_valley_lines = do_lines && do_valleys;
        do_ridge_surfaces = do_surfaces && do_ridges;
        do_valley_surfaces = do_surfaces && do_valleys;
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

void get_edge_reference(std::vector<edge_index_type>& edges, int cell_type) {
    if (cell_type == VTK_TETRA) {
        edges.push_back(edge_index_type(0, 1));
        edges.push_back(edge_index_type(1, 2));
        edges.push_back(edge_index_type(2, 0));
        edges.push_back(edge_index_type(0, 3));
        edges.push_back(edge_index_type(1, 3));
        edges.push_back(edge_index_type(2, 3));
    }
    else if (cell_type == VTK_HEXAHEDRON || cell_type == VTK_VOXEL) {
        edges.push_back(edge_index_type(0, 1));
        edges.push_back(edge_index_type(1, 2));
        edges.push_back(edge_index_type(2, 3));
        edges.push_back(edge_index_type(3, 0));
        edges.push_back(edge_index_type(4, 5));
        edges.push_back(edge_index_type(5, 6));
        edges.push_back(edge_index_type(6, 7));
        edges.push_back(edge_index_type(7, 4));
        edges.push_back(edge_index_type(0, 4));
        edges.push_back(edge_index_type(1, 5));
        edges.push_back(edge_index_type(2, 6));
        edges.push_back(edge_index_type(3, 7));
    }
    else if (cell_type == VTK_WEDGE) {
        edges.push_back(edge_index_type(0, 1));
        edges.push_back(edge_index_type(1, 2));
        edges.push_back(edge_index_type(2, 0));
        edges.push_back(edge_index_type(3, 4));
        edges.push_back(edge_index_type(4, 5));
        edges.push_back(edge_index_type(5, 3));
        edges.push_back(edge_index_type(0, 3));
        edges.push_back(edge_index_type(1, 4));
        edges.push_back(edge_index_type(2, 5));
    }
    else if (cell_type == VTK_PYRAMID) {
        edges.push_back(edge_index_type(0, 1));
        edges.push_back(edge_index_type(1, 2));
        edges.push_back(edge_index_type(2, 3));
        edges.push_back(edge_index_type(3, 0));
        edges.push_back(edge_index_type(0, 4));
        edges.push_back(edge_index_type(1, 4));
        edges.push_back(edge_index_type(2, 4));
        edges.push_back(edge_index_type(3, 4));
    }
    else {
        std::cerr << "Unsupported cell type: " << cell_type << '\n';
        exit(-1);
    }
}

void get_face_reference(std::vector<triangle_index_type>& triangles,
                        std::vector<quad_index_type>& quads,
                        int cell_type) {
    triangles.clear();
    quads.clear();
    if (cell_type == VTK_TETRA) {
        triangles.push_back(triangle_index_type(0, 1, 2));
        triangles.push_back(triangle_index_type(0, 1, 3));
        triangles.push_back(triangle_index_type(0, 2, 3));
        triangles.push_back(triangle_index_type(1, 2, 3));
    }
    else if (cell_type == VTK_HEXAHEDRON or cell_type == VTK_VOXEL) {
        quads.push_back(quad_index_type(0, 1, 2, 3));
        quads.push_back(quad_index_type(4, 5, 6, 7));
        quads.push_back(quad_index_type(0, 1, 5, 4));
        quads.push_back(quad_index_type(1, 2, 6, 5));
        quads.push_back(quad_index_type(2, 3, 7, 6));
        quads.push_back(quad_index_type(3, 0, 4, 7));
    }
    else if (cell_type == VTK_WEDGE) {
        triangles.push_back(triangle_index_type(0, 1, 2));
        triangles.push_back(triangle_index_type(3, 4, 5));
        quads.push_back(quad_index_type(0, 3, 4, 1));
        quads.push_back(quad_index_type(1, 4, 5, 2));
        quads.push_back(quad_index_type(2, 5, 3, 0));
    }
    else if (cell_type == VTK_PYRAMID) {
        quads.push_back(quad_index_type(0, 1, 2, 3));
        triangles.push_back(triangle_index_type(0, 1, 4));
        triangles.push_back(triangle_index_type(1, 2, 4));
        triangles.push_back(triangle_index_type(2, 3, 4));
        triangles.push_back(triangle_index_type(3, 0, 4));
    }
    else {
        std::cerr << "Unsupported cell type: " << cell_type << '\n';
        exit(-1);
    }
}

bool check_edge_linear(position_type& p, const edge_index_type& edge_id,
                       VTK_SMART(vtkDataSet) dataset,
                       const std::string& determinant_name) {
    // first time encountering this edge
    double value0, value1;
    VTK_SMART(vtkDataArray) determinant = dataset->GetPointData()->GetArray(determinant_name.c_str());
    determinant->GetTuple(edge_id[0], &value0);
    determinant->GetTuple(edge_id[1], &value1);
    if (value0 * value1 < 0) {
        // zero crossing found. compute linear approximation
        double u = -value0/(value1-value0);
        position_type p0, p1;
        dataset->GetPoint(edge_id[0], &p0[0]);
        dataset->GetPoint(edge_id[1], &p1[0]);
        p = (1.-u)*p0 + u*p1;
        return true;
    }
    else return false;
}

value_type PS_determinant(const vector_type& g, const matrix_type& H) {
    matrix_type A;
    A.column(0) = g;
    A.column(1) = H*g;
    A.column(2) = H*A.column(1);
    return determinant(A);
}

struct CreasePoint {
    CreasePoint(const position_type& p, const edge_index_type& eid,
                const vector_type& g, const matrix_type& H)
        : _p(p), _edge_id(eid), _g(g), _H(H) {
            
        spurt::sym_eigensystem(_evals, _evecs, _H);
        // crease type:
        //
    }
    int _crease_kind; // 0: ridge, 1: valley, -1: invalid
    int _crease_dim; // 1 or 2 dimension
    value_type _strength;
    position_type _p;
    edge_index_type _edge_id;
    vector_type _g;
    matrix_type _H;
    vec3 _evals;
    matrix_type _evecs;
};

struct CreaseInfo {
    CreaseInfo() : kind(-1), dimension(0), strength(0), position() {}
    int kind;
    int dimension;
    value_type strength;
    position_type position;
    value_type value;
};

CreaseInfo eberly_filter(const vector_type& g,
                         const matrix_type& evecs,
                         const vector_type& evals,
                         const double epsilon = eps) {
    // evals are sorted in increasing order
    vector_type g_dot_ev = abs(transpose(evecs)*g)/norm(g);
    // ridge test:
    int ridge_k=0; // ridge dimension
    for (int i=0; i<3; ++i, ++ridge_k) {
        if (g_dot_ev[i] > epsilon) break;
        if (evals[i] >= 0) break;
    }

    CreaseInfo info;
    if (ridge_k > 0) {
        info.kind = 0;
        info.dimension = ridge_k;
        info.strength = -evals[ridge_k-1];
        return info;
    }

    // valley test:
    int valley_k=0; // valley dimension
    for (int i=0; i<3; ++i, ++valley_k) {
        if (g_dot_ev[3-i] > epsilon) break;
        if (evals[i] <= 0) break;
    }
    if (valley_k > 0) {
        info.kind = 1;
        info.dimension = valley_k;
        info.strength = evals[valley_k-1];
        return info;
    }
    else return info;
}

bool linear_parallel_operator(std::vector<CreaseInfo>& cps,
                              const std::vector<vector_type>& g,
                              const std::vector<matrix_type>& H) {
    cps.clear();
    // finds locations where the linear vector fields g and H*g are parallel
    matrix_type V, W;
    /*
        v(s,t) = V(s, t, 1)^t
        v(b0, b1, b2) = b0 v0 + b1 v1 + b2 v2 = v0 + b1(v1-v0) + b2(v2-v0)
        V = [ v1x-v0x  v2x-v0x  v0x ][ b1 ]
            [ v1y-v0y  v2y-v0y  v0y ][ b2 ]
            [ v1z-v0z  v2z-v0z  v0z ][ 1  ]
    */

    std::vector<vector_type> Hg(3);
    for (unsigned int i = 0 ; i < 3 ; i++) {
        Hg[i] = H[i]*g[i];
    }
    V.column(0) = Hg[1]-Hg[0];
    V.column(1) = Hg[2]-Hg[0];
    V.column(2) = Hg[0];

    W.column(0) = g[1]-g[0];
    W.column(1) = g[2]-g[0];
    W.column(2) = g[0];

    double detW = std::abs(determinant(W));
    double detV = std::abs(determinant(V));
    matrix_type M;
    if (detW >= detV) {
        if (detW == 0) return false;
        M = inverse(W) * V;
    }
    else { // (detV > detW)
        M = inverse(V) * W;
    }
    
    complex_vector_type evals;
    complex_matrix_type evecs;
    eigensystem(evals, evecs, M);
    std::vector<position_type> found;
    for (int i=0; i<3; ++i) {
        if (evals[i].imag() != 0) break; // discard complex eigenvalues
        auto ev = evecs.column(i);
        value_type ez = ev[2].real();
        if (ez == 0.) continue; // degenerate solution
        ev /= ez; // divide by ev[2] to obtain homogeneous coordinates
        value_type b[3];
        b[1] = ev[0].real();
        b[2] = ev[1].real();
        b[0] = 1.-b[1]-b[2];
        if (b[1] >= 0 && b[2] >= 0 && b[0] >= 0) {
            // position lies inside the triangle
            found.push_back(position_type(b[0], b[1], b[2]));
        }
    }
    if (found.empty()) return false;
    else if (found.size() > 1) {
        std::cout << "Multiple PVO solutions found\n";
    }
    // loop over valid solutions inside triangles
    for (int i=0; i<found.size(); ++i) {
        position_type beta = found[i];
        matrix_type hX = beta[0]*H[0] + beta[1]*H[1] + beta[2]*H[2];
        vector_type gX = beta[0]*g[0] + beta[1]*g[1] + beta[2]*g[2];
        vector_type evalsX;
        matrix_type evecsX;
        sym_eigensystem(evalsX, evecsX, hX);
        CreaseInfo info = eberly_filter(gX, evecsX, evalsX, 0.01);
        if (info.kind != -1) {
            info.position = beta;
            cps.push_back(info);
            return true;
        }
    }
    return false;
}

// Interpolator
template<typename Interpolator>
struct cross_product_evaluator {
    typedef Interpolator interpolator_type;
    typedef typename interpolator_type::tensor_type tensor_type;

    cross_product_evaluator(const interpolator_type& intp)
        : m_interpolator(intp) {}

    bool value(vector_type& f, const position_type& p) const {
        vector_type g;
        matrix_type H;
        if (m_interpolator.gradient(g, p) &&
            m_interpolator.hessian(H, p)) {
            f = cross(g, H*g);
            return true;
        }
        return false;
    }
    bool gradient(matrix_type& df, const position_type& p) const {
        // derivative of f = g x Hg
        // f and g: vector fields, H: 2nd order tensor field
        // df = dg x Hg + g x ((dH)g + Hdg)
        // Note: dg/dxi = H.row(i)
        // Note: dH/dxi := dH[i]
        // df/dxi = H.row(i).cross(H*g) +
        //          g.cross(dH[i]*g + H*H.row(i))
        // needs g, H, dH
        vector_type g;
        matrix_type H;
        tensor_type dH;
        if (m_interpolator.gradient(g, p) &&
            m_interpolator.hessian(H, p) &&
            m_interpolator.dhessian(dH, p)) {
            for (int i=0; i<3; ++i) {
                df.row(i) = cross(H.row(i), H*g) +
                            cross(g, dH[i]*g + H*H.row(i));
            }
            return true;
        }
        return false;
    }

    interpolator_type m_interpolator;
};

//////////////////////////
// default interpolator for all non-BARG meshes
template<typename DataSet>
struct C3_vector_interpolator {
    typedef DataSet  dataset_type;
    typedef C3_vector_interpolator<dataset_type> self_type;
    typedef Eigen::VectorXf tuple_type;
    typedef typename vtk_utils::interpolator<DataSet, double, 3, vector_type, matrix_type, tuple_type>  sub_interpolator_type;
    typedef std::array<matrix_type, 3> tensor_type;

    C3_vector_interpolator(VTK_SMART(dataset_type) dataset,
                 const std::string& vector_name)
        : m_interpolator(dataset), m_dataset(dataset) {
        dataset->GetPointData()->SetActiveVectors(vector_name.c_str());
    }

    VTK_SMART(dataset_type) get_dataset() const {
        return m_dataset;
    }

    void set_max_order(int) {}

    bool value(vector_type& v, const position_type& x) const {
        if (m_interpolator.interpolate(v, x)) {
            return true;
        }
        return false;
    }

    bool value(vector_type& v, size_t index) const {
        m_dataset->GetPointData()->GetVectors()->GetTuple(index, &v[0]);
    }

    bool jacobian(matrix_type& J, const position_type& x) const {
        if (m_interpolator.interpolate(J, x)) {
            return true;
        }
        return false;
    }

    bool jacobian(matrix_type& J, size_t index) const {
        m_dataset->GetPointData()->GetArray("Jacobian")->GetTuple(index, &J(0,0));
        return true;
    }

    sub_interpolator_type m_interpolator;
    VTK_SMART(dataset_type) m_dataset;
};

template<>
struct C3_vector_interpolator<boundaryAwareRectGrid> {
    typedef boundaryAwareRectGrid dataset_type;
    typedef std::array<matrix_type, 3> tensor_type;
    typedef C3_vector_interpolator<dataset_type> self_type;

    C3_vector_interpolator(VTK_SMART(dataset_type) dataset,
                 const std::string& vector_name)
        : m_dataset(dataset), m_last_pos(invalid_value, 0, 0),
          m_max_order(3) {}

    VTK_SMART(dataset_type) get_dataset() const {
        return m_dataset;
    }

    void set_max_order(int order) {
        if (order != m_max_order && order > 0 && order < 4) {
            m_max_order = order;
            m_last_pos[0] = invalid_value;
        }
    }

    bool value(vector_type& v, const position_type& p) const {
        if (_interpolate(p)) {
            v = m_value;
            return true;
        }
        return false;
    }

    bool jacobian(matrix_type& J, const position_type& p) const {
        if (_interpolate(p)) {
            J = m_Jacobian;
            return true;
        }
        return false;
    }

    bool jacobian(matrix_type& J, size_t index) const {
        position_type p;
        m_dataset->GetPoint(index, &p[0]);
        if (_interpolate(p)) {
            J = m_Jacobian;
        }
        return true;
    }

    bool _interpolate(const position_type& p) const {
        // 0: f, 1: fx, 2: fy, 3: fz,
        // 4: fxx, 5: fxy, 6: fxz, 7: fyy, 8: fyz, 9: fzz,
        // 10: fxxx, 11: fxxy, 12: fxxz, 13: fxyy, 14: fxyz, 15: fxzz, 
        // 16: fyyy, 17: fyyz, 18: fyzz, 19: fzzz
        // dH[0] = [ [fxxx, fxxy, fxxz], [fxxy, fxyy, fxyz], 
        //           [fxxz, fxyz, fxzz] ]
        // dH[1] = [ [fxxy, fxyy, fxyz], [fxyy, fyyy, fyyz], 
        //           [fxyz, fyyz, fyzz] ]
        // dH[2] = [ [fxxz, fxyz, fxzz], [fxyz, fyyz, fyzz], 
        //           [fxzz, fyzz. fzzz] ]
        std::vector<value_type> value;
        if (is_bad_value(m_last_pos[0]) || any(p != m_last_pos)) {
            if (m_dataset->BsplineAllDerivatives(const_cast<double*>(&p[0]), value, 1, true)>=0) {
                // std::cout << "returned value:\n";
                // std::copy(value.begin(), value.end(), std::ostream_iterator<double>(std::cout, ","));
                // std::cout << '\n';
                // g = [fx, fy, fz]
                std::copy(&value[0], &value[3], &m_value[0]);
                std::copy(&value[3], &value[12], &m_Jacobian(0,0));
                return true;
            }
            else {
                // std::cerr << "Unable to interpolate at " << p << '\n';
                return false;
            }
            return false;
        }
        // else nothing to be done: derivatives at this position are
        // readily available
        return true;
    }

    mutable VTK_SMART(dataset_type) m_dataset;
    mutable position_type m_last_pos;
    mutable vector_type m_value;
    mutable matrix_type m_Jacobian;
    mutable int m_max_order;
};


//////////////////////////


// default interpolator for all non-BARG meshes
template<typename DataSet>
struct C3_interpolator {
    typedef DataSet  dataset_type;
    typedef C3_interpolator<dataset_type> self_type;
    typedef Eigen::VectorXf tuple_type;
    typedef typename vtk_utils::interpolator<DataSet, double, 3, vector_type, matrix_type, tuple_type>  sub_interpolator_type;
    typedef std::array<matrix_type, 3> tensor_type;

    C3_interpolator(VTK_SMART(dataset_type) dataset,
                 const std::string& scalar_name)
        : m_interpolator(dataset), m_dataset(dataset) {
        dataset->GetPointData()->SetActiveScalars(scalar_name.c_str());
    }

    VTK_SMART(dataset_type) get_dataset() const {
        return m_dataset;
    }

    void set_max_order(int) {}

    bool value(value_type& v, const position_type& x) const {
        if (m_interpolator.interpolate(v, x)) {
            return true;
        }
        return false;
    }

    bool gradient(vector_type& g, const position_type& x) const {
        if (m_interpolator.interpolate(g, x)) {
            return true;
        }
        return false;
    }

    bool gradient(vector_type& g, size_t index) const {
        m_dataset->GetPointData()->GetArray("Gradient")->GetTuple(index, &g[0]);
        return true;
    }

    bool hessian(matrix_type& H, const position_type& x) const {
        if (m_interpolator.interpolate(H, x)) {
            return true;
        }
        return false;
    }

    bool hessian(matrix_type& H, size_t index) const {
        m_dataset->GetPointData()->GetArray("Hessian")->GetTuple(index, &H(0,0));
        return true;
    }

    bool dhessian(tensor_type& dH, const position_type& x) const {
        Eigen::VectorXd value;
        if (m_interpolator.interpolate(value, x, "dHessian")) {
            // dH(0,0)/dx dH(0,0)/dy dH(0,0)/dz dH(1,0)/dx...
            for (int n=0; n<27; ++n) {
                int slice = n%3;
                int id = n/3;
                int i = id%3;
                int j = id/3;
                dH[slice](i,j) = value(n);
            }
            return true;
        }
        return false;
    }

    bool dHessian(tensor_type& dH, size_t index) const {
        double tuple[27];
        m_dataset->GetPointData()->GetArray("dHessian")->GetTuple(index, tuple);
        int counter=0;
        for (int col=0; col<3; ++col) {
            for (int row=0; row<3; ++row) {
                for (int slice=0; slice<3; ++slice) {
                    dH[slice](row,col) = tuple[counter++];
                }
            }
        }
        return true;
    }

    sub_interpolator_type m_interpolator;
    VTK_SMART(dataset_type) m_dataset;
};

template<>
struct C3_interpolator<boundaryAwareRectGrid> {
    typedef boundaryAwareRectGrid dataset_type;
    typedef std::array<matrix_type, 3> tensor_type;
    typedef C3_interpolator<dataset_type> self_type;

    C3_interpolator(VTK_SMART(dataset_type) dataset,
                 const std::string& scalar_name)
        : m_dataset(dataset), m_last_pos(invalid_value, 0, 0),
          m_max_order(3) {}

    VTK_SMART(dataset_type) get_dataset() const {
        return m_dataset;
    }

    void set_max_order(int order) {
        if (order != m_max_order && order > 0 && order < 4) {
            m_max_order = order;
            m_last_pos[0] = invalid_value;
        }
    }

    bool value(value_type& v, const position_type& p) const {
        if (_interpolate(p)) {
            v = m_value;
            return true;
        }
        return false;
    }

    bool gradient(vector_type& g, const position_type& p) const {
        if (_interpolate(p)) {
            g = m_gradient;
            return true;
        }
        return false;
    }

    bool gradient(vector_type& g, size_t index) const {
        position_type p;
        m_dataset->GetPoint(index, &p[0]);
        if (_interpolate(p)) {
            g = m_gradient;
        }
        return true;
    }

    bool hessian(matrix_type& H, const position_type& p) const {
        if (m_max_order>=2 && _interpolate(p)) {
            H = m_Hessian;
            return true;
        }
        return false;
    }

    bool hessian(matrix_type& H, size_t index) const {
        if (m_max_order < 2) return false;
        position_type p;
        m_dataset->GetPoint(index, &p[0]);
        if (_interpolate(p)) {
            H = m_Hessian;
        }
        return true;
    }

    bool dhessian(tensor_type& dH, const position_type& p) const {
        if (m_max_order>=3 && _interpolate(p)) {
            dH = m_dHessian;
            return true;
        }
        return false;
    }

    bool dHessian(tensor_type& dH, size_t index) const {
        if (m_max_order < 3) return false;
        position_type p;
        m_dataset->GetPoint(index, &p[0]);
        if (_interpolate(p)) {
            dH = m_dHessian;
        }
        return true;
    }

    bool _interpolate(const position_type& p) const {
        // 0: f, 1: fx, 2: fy, 3: fz,
        // 4: fxx, 5: fxy, 6: fxz, 7: fyy, 8: fyz, 9: fzz,
        // 10: fxxx, 11: fxxy, 12: fxxz, 13: fxyy, 14: fxyz, 15: fxzz, 16: fyyy, 17: fyyz, 18: fyzz, 19: fzzz
        // dH[0] = [ [fxxx, fxxy, fxxz], [fxxy, fxyy, fxyz], [fxxz, fxyz, fxzz] ]
        // dH[1] = [ [fxxy, fxyy, fxyz], [fxyy, fyyy, fyyz], [fxyz, fyyz, fyzz] ]
        // dH[2] = [ [fxxz, fxyz, fxzz], [fxyz, fyyz, fyzz], [fxzz, fyzz. fzzz] ]
        std::vector<value_type> value;
        if (is_bad_value(m_last_pos[0]) || any(p != m_last_pos)) {
            if (m_dataset->BsplineAllDerivatives(const_cast<double*>(&p[0]), value, m_max_order, true)>=0) {
                // std::cout << "returned value:\n";
                // std::copy(value.begin(), value.end(), std::ostream_iterator<double>(std::cout, ","));
                // std::cout << '\n';
                // g = [fx, fy, fz]
                m_value = value[0];
                std::copy(&value[1], &value[4], &m_gradient[0]);
                if (m_max_order >= 2) {
                    // H.row(0) = [fxx, fxy, fxz]
                    std::copy(&value[4], &value[7], &m_Hessian(0,0));
                    // symmetry: fxy = fyx
                    m_Hessian(0,1) = m_Hessian(1,0);
                    // symmetry: fxz = fzx
                    m_Hessian(0,2) = m_Hessian(2,0);
                    // H(1,1) = fyy
                    m_Hessian(1,1) = value[7];
                    // H(1,2) = fyz = fzy = H(2,1)
                    m_Hessian(1,2) = m_Hessian(2,1) = value[8];
                    // H(2,2) = fzz
                    m_Hessian(2,2) = value[9];
                    if (m_max_order > 2) {
                        // dH/dx(0,0) = fxxx
                        m_dHessian[0](0,0) = value[10];
                        // dH/dx(1,0) = fxyx = fxxy = dH/dy(0,0) = fyxx = dH/dx(0,1)
                        m_dHessian[0](1,0) = m_dHessian[0](0,1) = m_dHessian[1](0,0) = value[11];
                        // dH/dx(2,0) = fzxx = fxzx = dH/dx(0,2) = fxxz = dH/dz(0,0)
                        m_dHessian[0](2,0) = m_dHessian[0](0,2) = m_dHessian[2](0,0) = value[12];
                        // dH/dx(1,1) = fyyx = fxyy = dH/dy(0,1) = fyxy = dH/dy(1,0)
                        m_dHessian[0](1,1) = m_dHessian[1](0,1) = m_dHessian[1](1,0) = value[13];
                        // dH/dx(1,2) = fyzx = fzyx = dH/dx(2,1) = fzxy = dH/dy(2,0)
                        // = fxzy = dH/dy(0,2) = fxyz = dH/dz(0,1) = fyxz = dH/dz(1,0)
                        m_dHessian[0](1,2) = m_dHessian[0](2,1) = m_dHessian[1](0,2) = m_dHessian[1](2,0) = m_dHessian[2](0,1) = m_dHessian[2](1,0) = value[14];
                        // dH/dx(2,2) = fzzx = fzxz = dH/dz(2,0) = fxzz = dH/dz(0,2)
                        m_dHessian[0](2,2) = m_dHessian[2](2,0) = m_dHessian[2](0,2) = value[15];
                        // dH/dy(1,1) = fyyy
                        m_dHessian[1](1,1) = value[16];
                        // dH/dy(1,2) = fyzy = fzyy = dH/dy(2,1) = fyyz = dH/dz(1,1)
                        m_dHessian[1](1,2) = m_dHessian[1](2,1) = m_dHessian[2](1,1) = value[17];
                        // dH/dy(2,2) = fzzy = fzyz = dH/dz(2,1) = fyzz = dH/dz(1,2)
                        m_dHessian[1](2,2) = m_dHessian[2](2,1) = m_dHessian[2](1,2) = value[18];
                        // dH/dz(2,2) = fzzz
                        m_dHessian[2](2,2) = value[19];
                        m_last_pos = p;
                    } // order > 2
                } // order > 1
                // std::cout << "value=" << value[0] << '\n';
                // std::cout << "gradient=\n" << m_gradient << '\n';
                // std::cout << "hessian=\n" << m_Hessian << '\n';
                // std::cout << "dhessian/dx=\n" << m_dHessian[0] << '\n';
                // std::cout << "dhessian/dy=\n" << m_dHessian[1] << '\n';
                // std::cout << "dhessian/dz=\n" << m_dHessian[2] << '\n';
                return true;
            }
            else {
                // std::cerr << "Unable to interpolate at " << p << '\n';
                return false;
            }
            return false;
        }
        // else nothing to be done: derivatives at this position are
        // readily available
        return true;
    }

    mutable VTK_SMART(dataset_type) m_dataset;
    mutable position_type m_last_pos;
    mutable value_type m_value;
    mutable vector_type m_gradient;
    mutable matrix_type m_Hessian;
    mutable tensor_type m_dHessian;
    mutable int m_max_order;
};

template<typename DataSet>
using cross_product_type = cross_product_evaluator<C3_interpolator<DataSet>>;

template<typename Function>
bool lnsearch(const Function& func, position_type& x,
              vector_type& f0,
              const vector_type& dd, const value_type& maxlength)
{
    value_type lambda = 1.0;
    const value_type alpha = 1e-4;

    position_type xsave = x;
    vector_type fsave = f0;
    vector_type d = norm(dd) > maxlength ? dd * maxlength / norm(dd) : dd;
    vector_type f;

    for (unsigned int i = 0; i < 7; ++i) {
        x = xsave + lambda * d;
        if (func.value(f, x)) {
            if (norm(f) < (1 - alpha*lambda)*norm(fsave)) {
                return true;
            }
        }
        else {
            return false;
        }

        lambda *= 0.5;
    }

    return false;
}

template<typename Function>
bool newton_search(position_type& solution,
                   const std::vector<position_type>& face,
                   const Function& func, unsigned int nbiter=10,
                   value_type epsilon=eps) {
    // Newton algorithm with line search relaxation
    // f(x+dx) = 0
    // f(x+dx) = f(x) + df/dx(x) dx + O(dx^2) = 0
    // dx ~= -(df/dx)^{-1}*f(x)

    // intialize search at center
    position_type center = 0;
    for (int i=0; i<face.size(); ++i) {
        center += face[i];
    }
    center *= 1./value_type(face.size());
    position_type guess = center;
    value_type maxlength = 0.25*(norm(face[2]-face[0])+norm(face[3]-face[1]));

    // compute local 2d reference frame in cell
    vector_type e0, e1, normal;
    if (face.size() == 4) { // quadrilateral face
        e0 = face[2]-face[0];
        e1 = face[3]-face[1];
        normal = cross(e0, e1);
        e0 /= norm(e0);
        e1 = cross(normal,e0);
        e1 /= norm(e1);
    }

    // initialize search
    try {
        position_type x = guess;
        vector_type f;
        matrix_type df;
        if (!func.value(f, x) || !func.gradient(df, x)) {
            throw std::runtime_error("Unable to interpolate in Newton search");
        }

        for (int iter=0; iter<nbiter; ++iter) {
            matrix_type invdf = inverse(df);
            vector_type dx = -(invdf*f);
            // project dx on the linear approximation of the face
            dx = inner(dx, e0)*e0 + inner(dx, e1)*e1;
            if (!lnsearch(func, x, f, dx, maxlength)) {
                return false;
            }
            x += dx;
            if (!func.value(f, x) || !func.gradient(df, x)) {
                throw std::runtime_error("Unable to interpolate in Newton search");
            }
            if (norm(f) < epsilon) {
                solution = x;
                return true;
            }
        }
    }
    catch (std::exception& e) {
        std::cerr << "Newton search failed\n";
        std::cerr << e.what() << '\n';
    }
    return false;
}

template<typename Interpolator>
VTK_SMART(vtkDoubleArray) compute_determinant(Interpolator& intp) {
    std::cout << "computing determinant\n";
    size_t npoints = intp.get_dataset()->GetNumberOfPoints();
    VTK_CREATE(vtkDoubleArray, determinant);
    determinant->SetNumberOfComponents(1);
    determinant->SetNumberOfTuples(npoints);
    determinant->SetName("PS Determinant");
    intp.set_max_order(2);

    spurt::ProgressDisplay progress(true);
    progress.begin(npoints, "Computing determinant");
    for (size_t i=0; i<npoints; ++i) {
        vector_type g;
        matrix_type H;
        intp.gradient(g, i);
        intp.hessian(H, i);
        double det = PS_determinant(g, H);
        determinant->SetTypedTuple(i, &det);
        progress.update(i);
    }
    progress.end();
    return determinant;
}

template<typename Interpolator>
void compute_differential_quantities(Interpolator& intp,
    VTK_SMART(vtkDataSet) target, bool do_lambda2=true, bool do_vorticity=true, bool do_lamb=true, bool do_helicity=true) {
    VTK_SMART(vtkDataSet) dataset = target;
    size_t npoints = dataset->GetNumberOfPoints();
    VTK_CREATE(vtkDoubleArray, lambda2);
    VTK_CREATE(vtkDoubleArray, vorticity);
    VTK_CREATE(vtkDoubleArray, lamb);
    VTK_CREATE(vtkDoubleArray, helicity);
    if (do_lambda2) {
        lambda2->SetName("Lambda2");
        lambda2->SetNumberOfComponents(1);
        lambda2->SetNumberOfTuples(npoints);
    }
    if (do_lamb) {
        lamb->SetName("Lamb vector");
        lamb->SetNumberOfComponents(3);
        lamb->SetNumberOfTuples(npoints);
        do_vorticity = true;
    }
    if (do_helicity) {
        helicity->SetName("Helicity");
        helicity->SetNumberOfComponents(1);
        helicity->SetNumberOfTuples(npoints);
        do_vorticity = true;
    }
    if (do_vorticity) {
        vorticity->SetName("vorticity");
        vorticity->SetNumberOfComponents(3);
        vorticity->SetNumberOfTuples(npoints);
    }
    if (!do_lambda2 && !do_vorticity && !do_lamb && !do_helicity) {
        return;
    }
    for (size_t i=0; i<npoints; ++i) {
        vector_type velocity, w;
        matrix_type Jacobian;

        position_type p;
        dataset->GetPoint(i, &p[0]);
        if (intp.value(velocity, p) && intp.jacobian(Jacobian, p)) {
            if (do_vorticity) {
                w[0] = Jacobian(2,1) - Jacobian(1,2);
                w[1] = Jacobian(0,2) - Jacobian(2,0);
                w[2] = Jacobian(1,0) - Jacobian(0,1);
                vorticity->SetTuple3(i, w[0], w[1], w[2]);
            }
            if (do_lambda2) {
                matrix_type S = 0.5*(Jacobian + transpose(Jacobian));
                matrix_type Omega = 0.5*(Jacobian - transpose(Jacobian));
                S *= S;
                Omega *= Omega;
                S += Omega;
                vec3 evals;
                mat3 evecs;
                sym_eigensystem(evals, evecs, S);
                lambda2->SetTuple1(i, evals[1]);
            }
            if (do_lamb) {
                vector_type l = cross(velocity, w);
                lamb->SetTuple3(i, l[0], l[1], l[2]);
            }
            if (do_helicity) {
                value_type h = inner(velocity, w);
                helicity->SetTuple1(i, h);
            }
        }
    }
    if (do_vorticity)
        dataset->GetPointData()->AddArray(vorticity);
    if (do_lambda2)
        dataset->GetPointData()->AddArray(lambda2);
    if (do_lamb)
        dataset->GetPointData()->AddArray(lamb);
    if (do_helicity)
        dataset->GetPointData()->AddArray(helicity);
}

template<typename Interpolator>
void compute_crease_strength(Interpolator& intp) {
    std::cout << "computing crease strength\n";
    size_t npoints = intp.get_dataset()->GetNumberOfPoints();
    VTK_CREATE(vtkDoubleArray, rs);
    if (do_ridge_surfaces) {
        rs->SetNumberOfComponents(1);
        rs->SetNumberOfTuples(npoints);
        rs->SetName("Ridge surface strength");
    }
    VTK_CREATE(vtkDoubleArray, rl);
    if (do_ridge_lines) {
        rl->SetNumberOfComponents(1);
        rl->SetNumberOfTuples(npoints);
        rl->SetName("Ridge line strength");
    }
    VTK_CREATE(vtkDoubleArray, vs);
    if (do_valley_surfaces) {
        vs->SetNumberOfComponents(1);
        vs->SetNumberOfTuples(npoints);
        vs->SetName("Valley surface strength");
    }
    VTK_CREATE(vtkDoubleArray, vl);
    if (do_valley_lines) {
        vl->SetNumberOfComponents(1);
        vl->SetNumberOfTuples(npoints);
        vl->SetName("Valley line strength");
    }
    intp.set_max_order(2);

    spurt::ProgressDisplay progress(true);
    progress.begin(npoints, "Computing crease strength");
    for (size_t i=0; i<npoints; ++i) {
        // vector_type g;
        matrix_type H;
        // intp.gradient(g, i);
        vector_type eigenvalues;
        matrix_type eigenvectors;
        intp.hessian(H, i);
        sym_eigensystem(eigenvalues, eigenvectors, H);
        // eigenvalues are sorted in increasing order
        value_type ridge_surface_strength = std::max(-eigenvalues[2], 0.);
        value_type ridge_line_strength = std::max(-eigenvalues[1], 0.);
        value_type valley_surface_strength = std::max(eigenvalues[0], 0.);
        value_type valley_line_strength = std::max(eigenvalues[1], 0.);
        if (do_ridge_surfaces)
            rs->SetTypedTuple(i, &ridge_surface_strength);
        if (do_ridge_lines)
            rl->SetTypedTuple(i, &ridge_line_strength);
        if (do_valley_surfaces)
            vs->SetTypedTuple(i, &valley_surface_strength);
        if (do_valley_lines)
            vl->SetTypedTuple(i, &valley_line_strength);
        progress.update(i);
    }
    progress.end();
    if (do_ridge_surfaces)
        intp.get_dataset()->GetPointData()->AddArray(rs);
    if (do_ridge_lines)
        intp.get_dataset()->GetPointData()->AddArray(rl);
    if (do_valley_surfaces)
        intp.get_dataset()->GetPointData()->AddArray(vs);
    if (do_valley_lines)
        intp.get_dataset()->GetPointData()->AddArray(vl);
}

template<typename DataSet>
void filter_cells(std::vector<size_t>& valid_cells,
                  VTK_SMART(DataSet) dataset,
                  const std::string& scalar_name) {
    size_t ncells = dataset->GetNumberOfCells();
    valid_cells.clear();
    VTK_SMART(vtkDataArray) scalars = dataset->GetPointData()->GetArray(scalar_name.c_str());
    for (size_t i=0; i<ncells; ++i) {
        VTK_CREATE(vtkGenericCell, cell);
        dataset->GetCell(i, cell);
        int celltype = cell->GetCellType();
        VTK_SMART(vtkIdList) id_list = cell->GetPointIds();
        size_t nverts = id_list->GetNumberOfIds();
        bool valid = true;
        for (int k=0; k<nverts; ++k) {
            size_t vid = id_list->GetId(k);
            value_type v = scalars->GetTuple1(vid);
            value_type s;
            if (do_ridge_lines) {
                s = dataset->GetPointData()->GetArray("Ridge line strength")->GetTuple1(vid);
            }
            else if (do_ridge_surfaces) {
                s = dataset->GetPointData()->GetArray("Ridge surface strength")->GetTuple1(vid);
            }
            else if (do_valley_lines) {
                s = dataset->GetPointData()->GetArray("Valley line strength")->GetTuple1(vid);
            }
            else if (do_valley_surfaces) {
                s = dataset->GetPointData()->GetArray("Valley surface strength")->GetTuple1(vid);
            }
            if ((do_ridges && (v < value_threshold || s < strength_threshold)) ||
                (do_valleys && (v > value_threshold || s < strength_threshold))) {
                valid = false;
                break;
            }
        }
        if (valid) {
            valid_cells.push_back(i);
        }
    }
}

template<typename DataSet>
void extract_ridges(VTK_SMART(DataSet) dataset,
    const std::string& scalar_name) {
    typedef C3_interpolator<DataSet> interpolator_type;

    interpolator_type intp(dataset, scalar_name);

    bool changed = false;
    if (dataset->GetPointData()->GetArray("PS Determinant") == nullptr) {
        dataset->GetPointData()->AddArray(compute_determinant(intp));
        changed = true;
    }

    if ((do_ridge_surfaces &&
         dataset->GetPointData()->GetArray("Ridge surface strength") == nullptr)
        ||
        (do_ridge_lines &&
         dataset->GetPointData()->GetArray("Ridge line strength") == nullptr)
        ||
        (do_valley_surfaces &&
         dataset->GetPointData()->GetArray("Valley surface strength") == nullptr)
        ||
        (do_valley_lines &&
         dataset->GetPointData()->GetArray("Valley line strength") == nullptr)
        ) {
        compute_crease_strength(intp);
        changed = true;
    }

    std::vector<size_t> valid_cells;
    filter_cells(valid_cells, dataset, scalar_name);

    if (changed) {
        vtk_utils::saveVTK(dataset, name_in);
        std::cout << name_in << " dataset has been saved with the following fields:\n";
        int narrays = dataset->GetPointData()->GetNumberOfArrays();
        for (int i=0; i<narrays; ++i) {
            std::cout << "- " << dataset->GetPointData()->GetArray(i)->GetName() << '\n';
        }
    }

    // 3. Compute intersection points of 0-level sets of determinant
    std::map<edge_index_type, long int, lexicographical_order> edge_to_points;
    std::vector<position_type> all_crossings;
    std::vector<CreaseInfo> all_ridge_points;
    std::vector<CreaseInfo> all_valley_points;
    size_t ncells = valid_cells.size();
    std::cout << "There are " << ncells << " valid cells out of " << dataset->GetNumberOfCells() << " (" << 100.*ncells/dataset->GetNumberOfCells() << "%)\n";
    // First pass over all cells
    for (size_t i=0; i<ncells; ++i) {
        VTK_CREATE(vtkGenericCell, cell);
        dataset->GetCell(valid_cells[i], cell);
        int celltype = cell->GetCellType();
        std::vector<edge_index_type> edges;
        get_edge_reference(edges, celltype);

        VTK_SMART(vtkIdList) id_list = cell->GetPointIds();
        size_t nverts = id_list->GetNumberOfIds();
        // get global vertex indices
        std::vector<size_t> vert_ids;
        for (int k=0; k<nverts; ++k) {
            vert_ids.push_back(id_list->GetId(k));
        }
        // loop over cell edges
        for (int k=0; k<edges.size(); ++k) {
            auto edge = edges[k];
            edge_index_type edge_id(vert_ids[edge[0]], vert_ids[edge[1]]);
            if (edge_to_points.find(edge_id) == edge_to_points.end()) {
                position_type p;
                if (check_edge_linear(p, edge_id, dataset, "PS Determinant")) {
                        // NOT THREAD SAFE
                    all_crossings.push_back(p);
                    // reoord the index of the position we just computed
                    edge_to_points[edge_id] = all_crossings.size()-1;
                }
                else {
                    // no zero crossing
                    // assign invalid value to mark as checked
                    edge_to_points[edge_id] = -1;
                }
            }
        } // loop over edges
    } // First pass over all cells

    std::vector<CreaseInfo> crease_info(all_crossings.size());
    for (size_t i=0; i<crease_info.size(); ++i) {
        position_type p = all_crossings[i];
        vector_type g;
        matrix_type H;
        if (!intp.gradient(g, p) || !intp.hessian(H, p)) {
            std::cerr << "Unable to interpolate at " << p << "\n";
            continue;
        }
        vector_type eigenvalues;
        matrix_type eigenvectors;
        sym_eigensystem(eigenvalues, eigenvectors, H);
        crease_info[i] = eberly_filter(g, eigenvectors, eigenvalues);
    }

    size_t nactive_cells=0;
    for (size_t i=0; i<ncells; ++i) {
        VTK_CREATE(vtkGenericCell, cell);
        dataset->GetCell(valid_cells[i], cell);
        int celltype = cell->GetCellType();
        std::vector<edge_index_type> edges;
        std::vector<triangle_index_type> triangles;
        std::vector<quad_index_type> quads;
        get_edge_reference(edges, celltype);
        get_face_reference(triangles, quads, celltype);

        VTK_SMART(vtkIdList) id_list = cell->GetPointIds();
        size_t nverts = id_list->GetNumberOfIds();
        // get global vertex indices
        std::vector<size_t> vert_ids;
        for (int k=0; k<nverts; ++k) {
            vert_ids.push_back(id_list->GetId(k));
        }
        // loop over cell edges
/*
        std::vector<int> ridge_edges;
        std::vector<int> valley_edges;
        std::vector<position_type> ridge_points, valley_points;
        for (int k=0; k<edges.size(); ++k) {
            auto edge = edges[k];
            edge_index_type e(vert_ids[edge[0]], vert_ids[edge[1]]);
            if (edge_to_points[e] != -1) {
                const CreaseInfo& info = crease_info[edge_to_points[e]];
                if (info.kind == 0 && do_ridges) {
                    ridge_edges.push_back(k);
                    ridge_points.push_back(info.position);
                }
                else if (info.kind == 1 && do_valleys) {
                    valley_edges.push_back(k);
                    valley_points.push_back(info.position);
                }
            }
        }

        bool found_ridges = false;
        bool found_valleys = false;
        if (do_ridges && ridge_edges.size() >= 3) {
            std::cout << "Ridge surface detected in cell #" << valid_cells[i] << " (" << nverts << " vertices), " << ridge_edges.size() << " active edges\n";
            ++nactive_cells;
            found_ridges = true;
        }
        if (do_valleys && valley_edges.size() >= 3) {
            std::cout << "Valley surface detected in cell #" << valid_cells[i] << " (" << nverts << " vertices), " << valley_edges.size() << " active edges\n";
            ++nactive_cells;
            found_valleys = true;
        }
        if (ridge_edges.size()*valley_edges.size() > 0) {
            std::cout << "Ridge and valley types detected in same cell #" << valid_cells[i] << "(" << nverts << " vertices)!\n";
        }
*/

        if (do_lines) {
            if (true) { // || found_ridges || found_valleys) {
                // perform crease line point search on faces
                std::vector<CreaseInfo> ridge_points, valley_points;
                for (int t=0; t<triangles.size(); ++t) {
                    std::vector<vector_type> grads(3);
                    std::vector<matrix_type> hess(3);
                    std::vector<position_type> pts(3);
                    triangle_index_type tid = triangles[t];
                    for (int k=0; k<3; ++k) {
                        size_t id = vert_ids[tid[k]];
                        intp.gradient(grads[k], id);
                        intp.hessian(hess[k], id);
                        intp.get_dataset()->GetPoint(id, &pts[k][0]);
                    }
                    std::vector<CreaseInfo> found;
                    if (linear_parallel_operator(found, grads, hess)) {
                        std::cout << "Found crease point on triangle face #" << t << " in cell #" << valid_cells[i] << '\n';
                        // compute position
                        position_type bary = found[0].position;
                        position_type p = bary[0]*pts[0] + bary[1]*pts[1] + bary[2]*pts[2];
                        vector_type eigenvalues;
                        matrix_type eigenvectors;
                        matrix_type H;
                        intp.hessian(H, p);
                        vector_type g;
                        intp.gradient(g, p);
                        value_type v;
                        intp.value(v, p);
                        sym_eigensystem(eigenvalues, eigenvectors, H);
                        CreaseInfo info = eberly_filter(g, eigenvectors, eigenvalues);
                        if (do_ridges && info.kind==0) {
                            std::cout << "Found a ridge point with strength "
                                << -eigenvalues[1] << '\n';
                            info.position = p;
                            info.value = v;
                            info.strength = -eigenvalues[1];
                            ridge_points.push_back(info);
                        }
                        else if (do_valleys && info.kind == 1) {
                            std::cout << "Found a valley point with strength "
                                << eigenvalues[1] << '\n';
                            info.position = p;
                            info.value = v;
                            info.strength = eigenvalues[1];
                            valley_points.push_back(info);
                        }
                    }
                }
                if (valley_points.size() == 2) {
                    std::cout << "found a valley line in cell " << valid_cells[i]
                    << ": \n" << valley_points[0].position << "\n(" << valley_points[0].value << ")-(" << valley_points[1].value << ")\n" << valley_points[1].position << '\n';
                    all_valley_points.push_back(valley_points[0]);
                    all_valley_points.push_back(valley_points[1]);
                }
                else if (valley_points.size() > 0) {
                    std::cout << "valley points found but odd number: " << valley_points.size();
                }
                if (ridge_points.size() == 2) {
                    std::cout << "found a ridge line in cell " << valid_cells[i]
                    << ":\n" << ridge_points[0].position << "\n(" << ridge_points[0].value << ")-(" << ridge_points[1].value << ")\n" << ridge_points[1].position << '\n';
                    all_ridge_points.push_back(ridge_points[0]);
                    all_ridge_points.push_back(ridge_points[1]);
                }
                else if (ridge_points.size() > 0) {
                    std::cout << "ridge points found but odd number: " << ridge_points.size();
                }
            }
        }
    }
    // std::cout << "of the " << ncells << " valid cells, " << nactive_cells << " cells contain a crease (" << 100.*nactive_cells/ncells << "%)\n";

    if (all_valley_points.size() >= 2) {
        std::vector<position_type> points(all_valley_points.size());
        std::vector<value_type> values(all_valley_points.size());
        std::vector<value_type> strengths(all_valley_points.size());
        for (int i=0; i<all_valley_points.size(); ++i) {
            const CreaseInfo& info = all_valley_points[i];
            points[i] = info.position;
            values[i] = info.value;
            strengths[i] = info.strength;
        }

        VTK_SMART(vtkPolyData) lines = vtk_utils::make_points(points);
        vtk_utils::add_scalars(lines, values, true, "Value", false);
        vtk_utils::add_scalars(lines, strengths, true, "Valley strength", false);
        std::vector<size_t> indices(points.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::copy(indices.begin(), indices.end(), std::ostream_iterator<size_t>(std::cout, ", "));
        std::cout << '\n';
        vtk_utils::add_lines(lines, indices);
        vtk_utils::saveVTK(lines, name_out);
    }
}

int main(int argc, char* argv[]) {
    initialize(argc, argv);

    bool use_determinant = true; //do_surfaces || (do_lines && !use_newton);

    VTK_SMART(vtkDataSet) dataset = vtk_utils::readVTK(name_in);
    size_t npoints = dataset->GetNumberOfPoints();
    std::cout << "dataset imported. " << npoints << " points\n";

    // two cases: BARG dataset or other
    bool is_BARG = false;
    VTK_SMART(boundaryAwareRectGrid) BARG_dataset;
    if (vtkRectilinearGrid::SafeDownCast(dataset) != nullptr) {
        VTK_SMART(vtkRectilinearGrid) rgrid = vtkRectilinearGrid::SafeDownCast(dataset);
        // check if this is in fact a BARG dataset
        BARG_dataset = boundaryAwareRectGrid::New();
        BARG_dataset->ShallowCopy(rgrid);
        if (BARG_dataset->get_degree() >= 0) {
            is_BARG = true;
            BARG_dataset->SetArray();
            std::cout << "Dataset is BARG\n";
            std::cout << "Bspline degree is " << BARG_dataset->get_degree() << '\n';
        }
    }
    if (!is_BARG) {
        std::cout << "dataset is not BARG\n";
    }


#if 1 || DO_SCALARS
    std::cout << "scalar case selected at compile time\n";
    VTK_SMART(vtkDataArray) scalars;
    VTK_SMART(vtkDataArray) gradient;
    VTK_SMART(vtkDataArray) hessian;
    VTK_SMART(vtkDataArray) dhessian;

    scalars = dataset->GetPointData()->GetScalars();
    std::cout << "scalars = " << scalars << '\n';
    if (scalars == nullptr) {
        std::cerr << "ERROR: the dataset does not contain a scalar field\n";
        exit(-1);
    }
    std::cout << "scalar field has name " << scalars->GetName() << '\n';
    std::string scalar_name = scalars->GetName();

    if (!is_BARG) {
        bool changed_dataset = false;

        // 1. Compute / Import various derivatives we need
        gradient = dataset->GetPointData()->GetArray("Gradient");
        if (gradient.Get() == nullptr) {
            std::cout << "No gradient field present. Adding it\n";
            dataset = vtk_utils::add_gradient(dataset);
            std::cout << "Gradient field added\n";
            gradient = dataset->GetPointData()->GetArray("Gradient");
            changed_dataset = true;
        }
        else {
            std::cout << "Gradient field found\n";
        }
        hessian = dataset->GetPointData()->GetArray("Hessian");
        if (hessian.Get() == nullptr) {
            std::cout << "No Hessian field present. Adding it\n";
            dataset = vtk_utils::add_hessian(dataset);
            std::cout << "Hessian field added\n";
            hessian = dataset->GetPointData()->GetArray("Hessian");
            changed_dataset = true;
        }
        else {
            std::cout << "Hessian field found\n";
        }
        if (use_newton && do_lines) {
            dhessian = dataset->GetPointData()->GetArray("dHessian");
            if (dhessian.Get() == nullptr) {
                std::cout << "No Hessian derivative field present. Adding it\n";
                dataset = vtk_utils::add_derivative(dataset, "Hessian", "dHessian");
                std::cout << "Hessian derivative field added\n";
                dhessian = dataset->GetPointData()->GetArray("dHessian");
                changed_dataset = true;
            }
            else {
                std::cout << "Derivative of Hessian field found\n";
            }
        }
        if (changed_dataset) {
            vtk_utils::saveVTK(dataset, name_in);
            std::cout << "just exported " << name_in << " with 1st, 2nd, and 3rd derivatives added\n";
        }

        dataset->GetPointData()->SetActiveVectors("Gradient");
        dataset->GetPointData()->SetActiveTensors("Hessian");

        if (vtkUnstructuredGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkUnstructuredGrid) typed_data(vtkUnstructuredGrid::SafeDownCast(dataset));
            extract_ridges<vtkUnstructuredGrid>(typed_data, scalar_name);
        }
        else if (vtkRectilinearGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkRectilinearGrid) typed_data(vtkRectilinearGrid::SafeDownCast(dataset));
            extract_ridges<vtkRectilinearGrid>(typed_data, scalar_name);
        }
        else if (vtkStructuredGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkStructuredGrid) typed_data(vtkStructuredGrid::SafeDownCast(dataset));
            extract_ridges<vtkStructuredGrid>(typed_data, scalar_name);
        }
        else if (vtkUniformGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkUniformGrid) typed_data(vtkUniformGrid::SafeDownCast(dataset));
            extract_ridges<vtkUniformGrid>(typed_data, scalar_name);
        }
        else if (vtkImageData::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkImageData) typed_data(vtkImageData::SafeDownCast(dataset));
            extract_ridges<vtkImageData>(typed_data, scalar_name);
        }
        else {
            std::cerr << "Unsupported dataset type: " << dataset->GetClassName() << '\n';
            exit(-1);
        }
    }
    else {
        extract_ridges<boundaryAwareRectGrid>(BARG_dataset, scalar_name);
    }
#else
    std::cout << "Vector case selected at compile time\n"; 
    VTK_SMART(vtkDataArray) vectors;
    VTK_SMART(vtkDataArray) jacobian;

    VTK_SMART(vtkDataSet) surface;

    VTK_CREATE(vtkPlaneSource, plane);
    plane->SetOrigin(10000, -7500, 0);
    plane->SetXResolution(500);
    plane->SetYResolution(250);
    plane->SetPoint1(10000, 7500, 0);
    plane->SetPoint2(10000, -7500, 7500);
    plane->Update();
    surface = plane->GetOutput();
    surface->PrintSelf(std::cout, vtkIndent(0));

    // if (!name_surface.empty()) {
    //     surface = vtk_utils::readVTK(name_surface);
    // }
    // else {
    //     surface = dataset;
    // }

    vectors = dataset->GetPointData()->GetVectors();
    std::cout << "vectors = " << vectors << '\n';
    if (vectors == nullptr) {
        std::cerr << "ERROR: the dataset does not contain a vector field\n";
        exit(-1);
    }
    std::cout << "vector field has name " << vectors->GetName() << '\n';
    std::string vector_name = vectors->GetName();

    if (!is_BARG) {
        bool changed_dataset = false;

        // 1. Compute / Import various derivatives we need
        jacobian = dataset->GetPointData()->GetArray("Jacobian");
        if (jacobian.Get() == nullptr) {
            std::cout << "No Jacobian field present. Adding it\n";
            dataset = vtk_utils::add_jacobian(dataset);
            std::cout << "Jacobian field added\n";
            jacobian = dataset->GetPointData()->GetArray("Jacobian");
            changed_dataset = true;
        }
        else {
            std::cout << "Jacobian field found\n";
        }

        if (changed_dataset) {
            vtk_utils::saveVTK(dataset, name_in);
            std::cout << "just exported " << name_in << " with Jacobian added\n";
        }

        dataset->GetPointData()->SetActiveTensors("Jacobian");

        if (vtkUnstructuredGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkUnstructuredGrid) typed_data(vtkUnstructuredGrid::SafeDownCast(dataset));
            C3_vector_interpolator<vtkUnstructuredGrid> intp(typed_data, vectors->GetName());
            compute_differential_quantities(intp, surface, do_l2, do_omega, do_l, do_h);
        }
        else if (vtkRectilinearGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkRectilinearGrid) typed_data(vtkRectilinearGrid::SafeDownCast(dataset));
            C3_vector_interpolator<vtkRectilinearGrid> intp(typed_data, vectors->GetName());
            compute_differential_quantities(intp, surface, do_l2, do_omega, do_l, do_h);
        }
        else if (vtkStructuredGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkStructuredGrid) typed_data(vtkStructuredGrid::SafeDownCast(dataset));
            C3_vector_interpolator<vtkStructuredGrid> intp(typed_data, vectors->GetName());
            compute_differential_quantities(intp, surface, do_l2, do_omega, do_l, do_h);
        }
        else if (vtkUniformGrid::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkUniformGrid) typed_data(vtkUniformGrid::SafeDownCast(dataset));
            C3_vector_interpolator<vtkUniformGrid> intp(typed_data, vectors->GetName());
            compute_differential_quantities(intp, surface, do_l2, do_omega, do_l, do_h);
        }
        else if (vtkImageData::SafeDownCast(dataset) != nullptr) {
            VTK_SMART(vtkImageData) typed_data(vtkImageData::SafeDownCast(dataset));
            C3_vector_interpolator<vtkImageData> intp(typed_data, vectors->GetName());
            compute_differential_quantities(intp, surface, do_l2, do_omega, do_l, do_h);
        }
        else {
            std::cerr << "Unsupported dataset type: " << dataset->GetClassName() << '\n';
            exit(-1);
        }
    }
    else {
        C3_vector_interpolator<boundaryAwareRectGrid> intp(
            BARG_dataset, "velocity");
        compute_differential_quantities(intp, surface, do_l2, do_omega, do_l, do_h);
    }
    vtk_utils::saveVTK(surface, name_surface);
    std::cout << "Saved " << name_surface << " with following attributes included\n";
    for (int i=0; i<surface->GetPointData()->GetNumberOfArrays(); ++i) {
        std::cout << "- " << surface->GetPointData()->GetArray(i)->GetName() << '\n';
        surface->GetPointData()->GetArray(i)->PrintSelf(std::cout, vtkIndent(0));
    }
#endif

    return 0;
}
