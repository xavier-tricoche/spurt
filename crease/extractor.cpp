#include <cmath>
#include <crease/extractor.hpp>
#include <math/math.hpp>
#include <misc/sort.hpp>
#include <math/matrix.hpp>
#include <crease/grid.hpp>
#include <set>
#include <list>
#include <algorithm>
#include <fstream>
#include <teem/hest.h>
#include <misc/progress.hpp>
#include <crease/pvo.hpp>

// for debugging display
std::vector< std::vector< nvis::vec3 > > xavier::crease::vertices;
std::vector< nvis::vec3 > xavier::crease::problematic_voxels;
std::vector< nvis::vec3 > xavier::crease::fixed_voxels;
std::vector< nvis::vec3 > xavier::crease::current_vertices;
std::vector< nvis::vec3 > xavier::crease::ok_faces;
std::vector< unsigned int > xavier::crease::added_vertices;
std::vector< nvis::vec3 > xavier::crease::crossing_faces;
std::vector< xavier::crease::point_on_face > xavier::crease::all_points_on_face;
std::vector< xavier::crease::path > xavier::crease::paths;
std::vector< nvis::vec3 > xavier::crease::pvo_faces;
std::map< xavier::FaceId, nvis::vec3 > xavier::crease::reference_points;
xavier::FaceId xavier::crease::current_face_id;
std::vector< nvis::vec3 > xavier::crease::show_all;
std::vector< nvis::vec3 > xavier::crease::intermediate_steps;
std::vector< nvis::vec3 > xavier::crease::round1, xavier::crease::round2, xavier::crease::round12;
std::string xavier::crease::flag_file_name;
bool xavier::crease::bold_move;

// extraction control
double xavier::crease::value_threshold = 0;
double xavier::crease::strength_threshold = 0;
double xavier::crease::value_threshold_select = 0;
double xavier::crease::strength_threshold_select = 0;
double xavier::crease::confidence_threshold = 0.5;
double xavier::crease::gradient_eps = 1.0e-6;
double xavier::crease::gradient_eps_rel = 0.05;
bool xavier::crease::apply_filter;
bool xavier::crease::speedup;
bool xavier::crease::fixing_voxel;
bool xavier::crease::display_debug_info = false;
bool xavier::crease::read_info;

// spatial accuracy
unsigned int xavier::crease::max_depth = 4;
unsigned int xavier::crease::max_depth_fix = 6;
unsigned int xavier::crease::upsample;
double xavier::crease::max_int_error = 0.05;
double xavier::crease::max_align_error = 0.01;

unsigned int xavier::crease::nb_crossings;

// crease type
bool xavier::crease::is_ridge;
int xavier::crease::crease_kind;
unsigned int xavier::crease::failed_conv;
xavier::crease::extraction_method xavier::crease::ext_meth;

// interface to gage
xavier::MeasureWrapper* xavier::crease::the_wrapper;

using namespace xavier;

// global face information
xavier::crease::face_information xavier::crease::current_face_info;

namespace voxel_info {
// multipurpose voxel information
// faces
const unsigned int faces[6][4] = {
    { 0, 1, 2, 3 },
    { 4, 7, 6, 5 },
    { 0, 4, 5, 1 },
    { 1, 5, 6, 2 },
    { 2, 6, 7, 3 },
    { 3, 7, 4, 0 }
};

// procedural vertices
const unsigned int idx[8][3] = {
    { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
    { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 }
};

// procedural face vertices
unsigned int face_pt[3][4][3] = {
    { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 } },
    { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 } },
    { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } }
};

};

std::ostream& xavier::crease::operator<<(std::ostream& out, const face_type& face)
{
    out << "face: " << std::endl
    << "positions: 0:" << face.p[0] << ", 1:" << face.p[1]
    << ", 2:" << face.p[2] << ", 3:" << face.p[3] << std::endl;
    /*
        << "gradients: 0:" << face.g[0] << ", 1:" << face.g[1]
            << ", 2:" << face.g[2] << ", 3:" << face.g[3] << std::endl
            << "eigenvecs: 0:" << face.e[0] << ", 2:" << face.e[1]
            << ", 2:" << face.e[2] << ", 3:" << face.e[3] << std::endl;
    */
    return out;
}

unsigned int write_vertex_info(const Grid& grid, std::vector< bool >& ok)
{
    using namespace xavier;
    using namespace xavier::crease;

    std::ofstream output(xavier::crease::flag_file_name.c_str(), std::ios::binary);
    unsigned int nb_ok_vertices = 0;

    unsigned int M, N, P;
    M = grid.size[0];
    N = grid.size[1];
    P = grid.size[2];

    ok.resize(M*N*P);

    int prev_pct = -1;

    double *buffer = (double *)calloc(2 * M * N * P, sizeof(double));

    for (unsigned int k = 0 ; k < P ; ++k) {
        for (unsigned int j = 0 ; j < N ; ++j) {
            for (unsigned int i = 0 ; i < M ; ++i) {

                unsigned int id = i + M * (j + N * k);

                int pct = 100 * id / (M * N * P);
                if (pct > prev_pct) {
                    std::cout << "progress: " << pct << "% complete" << std::endl;
                    prev_pct = pct;
                }

                nvis::vec3 q = grid(i, j, k);
                double val = the_wrapper->value(q);
                double str = the_wrapper->eigenvalue(q, 1);

                ok[id] = (crease::is_ridge ? (val > value_threshold && str < strength_threshold) :
                          (val < value_threshold && str > strength_threshold));

                if (ok[id]) ++nb_ok_vertices;

                buffer[2*id] = val;
                buffer[2*id+1] = str;
            }
        }
    }

    Nrrd *nrrd = nrrdNew();
    if (nrrdWrap_va(nrrd, buffer, nrrdTypeDouble, 4, 2, M, N, P) ||
        nrrdSave(xavier::crease::flag_file_name.c_str(), nrrd, NULL)) {
        char *err = biffGetDone(NRRD);
        std::cout << err << std::endl;
    }

    nrrdNix(nrrd);
    delete[] buffer;

    return nb_ok_vertices;
}

int compute_vertex_info(const Grid& grid, std::vector< bool >& ok, std::vector< bool >& conf_ok)
{
    using namespace xavier;
    using namespace xavier::crease;

    std::cout << "computing vertex info..." << std::endl;
    std::cout << "confidence threshold = " << xavier::crease::confidence_threshold << std::endl;
    std::cout << "value threshold = " << value_threshold << std::endl;
    std::cout << "strength threshold = " << strength_threshold << std::endl;

    unsigned int M, N, P;
    M = grid.size[0];
    N = grid.size[1];
    P = grid.size[2];

    ok.resize(M*N*P);
    conf_ok.resize(M*N*P);

    unsigned int nb_ok_vertices = 0;

    nvis::vec3 p;
    double val, conf, str;

    int prev_pct = -1;
    for (unsigned int k = 0 ; k < P ; ++k) {
        for (unsigned int j = 0 ; j < N ; ++j) {
            for (unsigned int i = 0 ; i < M ; ++i) {

                unsigned int id = i + M * (j + N * k);

                int pct = 100 * id / (M * N * P);
                if (pct > prev_pct) {
                    std::cout << "progress: " << pct << "% complete" << std::endl;
                    prev_pct = pct;
                }

                p = grid(i, j, k);
                conf = the_wrapper->confidence(p);
                val = the_wrapper->value(p);
                str = the_wrapper->eigenvalue(p, 1);
                // std::cout << "(" << i << ", " << j << ", " << k << "): " << conf << "; " << val
                // << "; " << str << std::endl;

                conf_ok[id] = (conf > xavier::crease::confidence_threshold);
                // if (!conf_ok[id]) {
                //     std::cout << "low confidence at (" << i << ", " << j << ", " << k << ")" << std::endl;
                // }
                ok[id] =  conf_ok[id] &&
                          (crease::is_ridge ?
                           (val > value_threshold && str < strength_threshold) :
                           (val < value_threshold && str > strength_threshold));

                // if (ok[id]) {
                //     std::cout << "(" << i << ", " << j << ", " << k << ") is valid: "
                //     << conf << ", " << val << ", " << str << std::endl;
                // }

                if (ok[id]) ++nb_ok_vertices;
            }
        }
    }
    return nb_ok_vertices;
}

unsigned int read_vertex_info(std::vector< bool >& ok, std::vector< bool >& conf_ok)
{
    using namespace xavier;
    using namespace xavier::crease;

    Nrrd *nrrd = nrrdNew();
    if (nrrdLoad(nrrd, xavier::crease::flag_file_name.c_str(), NULL)) {
        char *err = biffGetDone(NRRD);
        std::cout << err << std::endl;
    }

    xavier::nrrd_utils::nrrd_data_wrapper<double> buffer(nrrd);

    bool has_confidence = (nrrd->dim == 4 && nrrd->axis[0].size == 3);

    unsigned int M, N, P;
    M = nrrd->axis[1].size;
    N = nrrd->axis[2].size;
    P = nrrd->axis[3].size;

    ok.resize(M*N*P);

    unsigned int nb_ok_vertices = 0;

    int prev_pct = -1;
    for (unsigned int k = 0 ; k < P ; ++k) {
        for (unsigned int j = 0 ; j < N ; ++j) {
            for (unsigned int i = 0 ; i < M ; ++i) {

                unsigned int id = i + M * (j + N * k);

                int pct = 100 * id / (M * N * P);
                if (pct > prev_pct) {
                    std::cout << "progress: " << pct << "% complete" << std::endl;
                    prev_pct = pct;
                }

                double val, str, cfd;
                if (has_confidence) {
                    val = buffer[3*id+1];
                    str = buffer[3*id+2];
                    cfd = buffer[3*id];
                }
                else {
                    val = buffer[2*id];
                    str = buffer[2*id+1];
                }

                conf_ok[id] = true;
                if (has_confidence) {
                    conf_ok[id] = (cfd > xavier::crease::confidence_threshold);
                }
                ok[id] =  conf_ok[id] &&
                          (crease::is_ridge ?
                           (val > value_threshold && str < strength_threshold) :
                           (val < value_threshold && str > strength_threshold));

                if (ok[id]) ++nb_ok_vertices;
            }
        }
    }
    nrrdNuke(nrrd);

    return nb_ok_vertices;
}

// ---------------------------------------------------------------------------

// compute PCA of a set of non-oriented vectors
void crease::PCA(std::vector< nvis::vec3 >& evec, const std::vector< nvis::vec3 >& dirs)
{
    double mat[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

// mat = sum_i( dirs[i] * dirs[i]^T )
    for (unsigned int l = 0 ; l < dirs.size() ; l++) {
        for (unsigned int i = 0 ; i < 3 ; i++)
            for (unsigned int j = i ; j < 3 ; j++)
                mat[i+3*j] += dirs[l][i] * dirs[l][j];
    }
    mat[3] = mat[1];
    mat[6] = mat[2];
    mat[7] = mat[5];

    double _eval[3];
    double _evec[9];
    int nbroots = ell_3m_eigensolve_d(_eval, _evec, mat, 1);
    evec.resize(3);
    for (unsigned int i = 0 ; i < 3 ; i++) {
        for (unsigned int j = 0 ; j < 3 ; j++) {
            evec[i][j] = _evec[3*i+j];
        }
    }
}

// Stetten's average direction
nvis::vec3 crease::average(const std::vector< nvis::vec3 >& dirs)
{
    std::vector< nvis::vec3 > tmp(3);
    PCA(tmp, dirs);
    return tmp[0];
}

// ---------------------------------------------------------------------------
// Furst's eigenvector orientation method
void crease::Furst_orientation(nvis::vec3& ev1, const nvis::vec3& ev0)
{
    std::vector< nvis::vec3 > vecs(2);
    vecs[0] = ev0;
    vecs[1] = ev1;
    nvis::vec3 avg = average(vecs);
    double dot0 = nvis::inner(ev0, avg);
    double dot1 = nvis::inner(ev1, avg);
    if (dot0*dot1 < 0)
        ev1 *= -1;
}

// ---------------------------------------------------------------------------
// Naive eigenvector orientation method
void crease::Naive_orientation(nvis::vec3& ev1, const nvis::vec3& ev0)
{
    if (nvis::inner(ev0, ev1) < 0)
        ev1 *= -1;
}

// ---------------------------------------------------------------------------

// 0: OK
// 1: invalid value(s)
// 2: invalid strength(s)
int crease::filter(const std::vector< double >& vals,
                   const std::vector< double >& strs)
{
    if (!apply_filter) return 0;

    assert(vals.size() == strs.size());

    bool gv = false; // is there one good value?
    bool bv = false; // is there one bad value?
    bool gs = false; // is there one valid strength?
    bool bs = false; // is there one invalid strength?

    for (unsigned int i = 0 ; i < vals.size() ; i++) {
        const double& val = vals[i];
        const double& strength = strs[i];

        bool val_ok = good_value(val);
// gv = (gv || val_ok);
        bv = (bs || !val_ok);

        bool str_ok = good_strength(strength);
        gs = (gs || str_ok);
// bs = (bs || !str_ok);
    }

// // for now conservative answer
// if (!gv) return 1;
// else if (!gs) return 2;
// else return 0;

// aggressive filtering: need all good values and all good strengths
    if (bv) return 1;
    else if (bs) return 2;
    else return 0;
}

// 0: OK
// 1: invalid value(s)
// 2: invalid strength(s)
// 3: invalid confidence
int crease::filter(const std::vector< double >& vals,
                   const std::vector< double >& strs,
                   const std::vector< double >& cfds)
{
    if (!apply_filter) return 0;

    assert(vals.size() == strs.size() && vals.size() == cfds.size());

    bool gv = false; // is there one good value?
    bool bv = false; // is there one bad value?
    bool gs = false; // is there one valid strength?
    bool bs = false; // is there one invalid strength?

    for (unsigned int i = 0 ; i < vals.size() ; i++) {
        if (cfds[i] < crease::confidence_threshold) return 3;

        const double& val = vals[i];
        const double& str = strs[i];

        bool val_ok = (crease::is_ridge ?
                       val > crease::value_threshold :
                       val < crease::value_threshold);
        gv = (gv || val_ok);
        bv = (bs || !val_ok);

        bool str_ok = (crease::is_ridge ?
                       str < crease::strength_threshold :
                       str > crease::strength_threshold);
        gs = (gs || str_ok);
        bs = (bs || !str_ok);
    }

// for now conservative answer
    if (!gv) return 1;
    else if (!gs) return 2;
    else return 0;
}

// --------------------------------------------------------------------------

inline FaceId face_id(unsigned int id, unsigned int axis,
                      unsigned int m, unsigned int n)
{
    unsigned int incr_x = 1;
    unsigned int incr_y = m;
    unsigned int incr_z = m * n;
    switch (axis) {
    case 0: {
// YZ
        return FaceId(id, id + incr_y, id + incr_y + incr_z, id + incr_z);
    }
    case 1: {
// XZ
        return FaceId(id, id + incr_x, id + incr_x + incr_z, id + incr_z);
    }
    case 2: {
// XY
        return FaceId(id, id + incr_x, id + incr_x + incr_y, id + incr_y);
    }
    default: { throw std::runtime_error("Invalid index in face_id()"); }
    }
}

inline void ijk(unsigned int& i, unsigned int& j, unsigned& k, unsigned int id,
                unsigned int M, unsigned int N)
{
    k = id / M / N;
    j = (id - k * M * N) / M;
    i = id % M;
}

// --------------------------------------------------------------------------

std::map< FaceId, unsigned int > face_found;
std::vector< nvis::vec3 > xavier::crease::all_face_points;
std::vector< std::pair< unsigned int, unsigned int > > all_edges;
unsigned int nb_segments;

// --------------------------------------------------------------------------

inline void add_segment(std::vector< unsigned int > ids)
{
    if (ids.size() < 2) return;
    else if (ids.size() == 2) {
        ++nb_segments;
        all_edges.push_back(std::pair< unsigned int, unsigned int >
                            (ids[0], ids[1]));
    }
    else if (ids.size() >= 3) {
        unsigned int np = ids.size();
        std::vector< double > _vals(np);
        for (unsigned int l = 0 ; l < np ; l++) {
            _vals[l] = crease::the_wrapper->value(crease::all_face_points[ids[l]]);
        }
        std::vector< unsigned int > _ids(np);
        xavier::sort_ids(_ids, _vals);
        ++nb_segments;
        unsigned int i0 = _ids[np-2];
        unsigned int i1 = _ids.back();
        all_edges.push_back(std::pair< unsigned int, unsigned int >
                            (ids[i0], ids[i1]));
    }
    std::cout << "found a segment" << std::endl;
}

nvis::vec3 trilinear(const nvis::vec3& p, const std::vector< nvis::vec3 >& v)
{
    nvis::vec3 u;
    const nvis::vec3& min = v[0];
    const nvis::vec3& max = v[6];
    if (xavier::crease::display_debug_info)
        std::cout << "min = " << min << ", max = " << max << std::endl;
    for (unsigned int i = 0; i < 3 ; ++i) {
        u[i] = (p[i] - min[i]) / (max[i] - min[i]);
    }
    return u;
}

bool inside(const nvis::vec3& local)
{
    return (local[0] > 0 && local[0] < 1 &&
            local[1] > 0 && local[1] < 1 &&
            local[2] > 0 && local[2] < 1);
}

void intersect(nvis::vec3& xp, unsigned int& fid,
               const nvis::vec3& p0, const nvis::vec3& p1)
{
    static const unsigned int table[][2] = { { 5, 3 }, { 2, 4 }, { 0, 1 } };

    double l = nvis::norm(p1 - p0);
    double t[] = { 0, 0, 0 };
    int c[] = { 0, 0, 0 };
    for (unsigned int i = 0 ; i < 3 ; i++) {
        if (p1[i] <= 0) {
            t[i] = -p0[i] / (p1[i] - p0[i]);
            c[i] = -1;
        }
        else if (p1[i] >= 1) {
            t[i] = (1 - p0[i]) / (p1[i] - p0[i]);
            c[i] = 1;
        }
    }
    double tmin = 2;
    int minid = -1;
    for (unsigned int i = 0 ; i < 3 ; i++) {
        if (t[i] != 0 && t[i] < tmin) {
            tmin = t[i];
            minid = i;
        }
    }
    xp = (1 - t[minid]) * p0 + t[minid] * p1;
    int k = (c[minid] < 0 ? 0 : 1);
    xp[minid] = k;
    fid = table[minid][k];
}

// --------------------------------------------------------------------------

bool track_ridge_line(nvis::vec3& out, unsigned int& fid_out,
                      const std::vector< nvis::vec3 >& voxel,
                      const nvis::vec3& entry,
                      const unsigned int fid_in)
{
    // get initial location
    nvis::vec3 start = entry;
    double step = 0.05 * nvis::norm(voxel[6] - voxel[0]); // 5% of voxel diagonal
    // initialize orientation
    unsigned int f = fid_in;
    const unsigned int* face = voxel_info::faces[f];
    nvis::vec3 fp[4];
    for (unsigned int l = 0 ; l < 4 ; l++) {
        fp[l] = voxel[face[l]];
    }
    // inward pointing face normal
    nvis::vec3 ref = nvis::cross(fp[1] - fp[0], fp[3] - fp[0]);

    // Euler integration along eigenvector
    nvis::vec3 p0 = start;
    nvis::vec3 l0 = trilinear(p0, voxel);
    nvis::vec3 p1;
    unsigned int n = 0;
    // my_path.push_back(p0);
    while (n < 100) {
        ++n;
        xavier::crease::intermediate_steps.push_back(p0);
        nvis::vec3 dir = xavier::crease::the_wrapper->eigenvector(p0, crease::is_ridge ? 0 : 2);
        if (nvis::inner(dir, ref) < 0) dir *= -1;
        p1 = p0 + step * dir; // Euler step in physical space
        nvis::vec3 l1 = trilinear(p1, voxel); // corresponding logical coordinates

        // did we leave the voxel?
        if (!inside(l1)) {
            nvis::vec3 stop;
            intersect(stop, fid_out, l0, l1);
            out = voxel[0] +
                  stop[0] * (voxel[1] - voxel[0]) +
                  stop[1] * (voxel[3] - voxel[0]) +
                  stop[2] * (voxel[4] - voxel[0]);
            return true;
        }
        ref = p1 - p0; // last step becomes reference direction
        p0 = p1; // step forward in physical space
        l0 = l1; // and in logical space
    }

    return false;
}

// --------------------------------------------------------------------------

bool trace_ray(nvis::vec3& out, unsigned int& fid_out,
               const std::vector< nvis::vec3 >& voxel,
               const nvis::vec3& entry,
               const nvis::vec3& ray,
               const unsigned int fid_in)
{
// determine step length that ensures crossing
    double diag = nvis::norm(voxel[0] - voxel[6]);
// solve for intersection point in local coordinates
    nvis::vec3 uvw0 = trilinear(entry, voxel);
    nvis::vec3 uvw1 = trilinear(entry + 1.5 * diag * ray, voxel);
    nvis::vec3 x;
    intersect(x, fid_out, uvw0, uvw1);
// convert to spatial coordinates
    out = voxel[0] +
          x[0] * (voxel[1] - voxel[0]) +
          x[1] * (voxel[3] - voxel[0]) +
          x[2] * (voxel[4] - voxel[0]);

    return true;
}

// --------------------------------------------------------------------------

FaceId global_face_id(unsigned int i, unsigned int j, unsigned int k,
                      const Grid& grid, unsigned int fid)
{
    const unsigned int * face = voxel_info::faces[fid];
    unsigned int id[4];
    for (unsigned int l = 0 ; l < 4 ; l++) {
        id[l] = grid.id(i + voxel_info::idx[face[l]][0],
                        j + voxel_info::idx[face[l]][1],
                        k + voxel_info::idx[face[l]][2]);
    }
    return FaceId(id[0], id[1], id[2], id[3]);
}

// --------------------------------------------------------------------------

bool check_confidence(unsigned int i, unsigned int j, unsigned int k,
                      unsigned int M, unsigned int N,
                      const std::vector< bool >& conf_ok)
{
    unsigned int x, y, z, id;
    for (unsigned int n = 0 ; n < 8 ; ++n) {
        x = i + voxel_info::idx[n][0];
        y = j + voxel_info::idx[n][1];
        z = k + voxel_info::idx[n][2];
        id = x + M * (y + N * z);
        if (!conf_ok[id]) return false;
    }
    return true;
}

// --------------------------------------------------------------------------

void crease::extract_lines(std::vector< line >& creases, const Nrrd* nrrd)
{
    using namespace voxel_info;

    struct pvo_solution {
        pvo_solution() : p(), local_fid(0), global_fid() {}

        nvis::vec3 p;
        unsigned int local_fid;
        FaceId global_fid;
    };


    std::cout << "we are extracting "
              << (is_ridge ? "RIDGES" : "VALLEYS") << std::endl;
    std::cout << "crease_kind = " << crease_kind << std::endl;
    std::cout << "spatial upsampling: " << crease::upsample << std::endl;
    std::cout << "subdividion depth = " << max_depth << std::endl;
    std::cout << "threshold on value = " << value_threshold << std::endl;
    std::cout << "threshold on crease strength = " << strength_threshold << std::endl;

    apply_filter = true;
    nb_crossings = 0;

    xavier::crease::search_face_PVO search_face;

    // check if this is a tensor field <-- do we need that?
    bool is_tensor = (nrrd->dim == 4 && nrrd->axis[0].size == 7);

    Grid sample_grid(nrrd, crease::upsample);
    std::cerr << "after initialization: grid bounding box is " << sample_grid.bounds
              << ", step = " << sample_grid.d << ", size " << sample_grid.size << '\n';

    std::set< FaceId > face_done; // all the voxel faces processed so far
    all_face_points.clear(); // all the crease points found
    all_edges.clear(); // all the crease segments found
    face_found.clear(); // all the faces where crease points were found so far
    ok_faces.clear();
    added_vertices.clear();

    failed_conv = 0;

    unsigned int nbok = 0;
    nb_segments = 0;
    unsigned int nb_several = 0;

    unsigned int M, N, P;
    M = sample_grid.size[0];
    N = sample_grid.size[1];
    P = sample_grid.size[2];

    unsigned int nb_wrong_val = 0;
    unsigned int nb_wrong_eval = 0;

    // voxel data
    std::vector< nvis::vec3 > v(4); // vertex coordinates
    std::vector< double > val(4);   // associated values
    std::vector< double > str(4);   // associated strength
    std::vector< double > cfd(4);   // associated confidence
    std::vector< unsigned int > ids(4);  // vertex indices

    crease::vertices.clear();

    const unsigned int invalid_face_id = std::numeric_limits< unsigned int >::max();
    std::map< FaceId, unsigned int > face_point_ids;

    // create a list of all the faces that need to be processed
    unsigned int incr_x = 1;
    unsigned int incr_y = M;
    unsigned int incr_z = M * N;

    unsigned int n1 = P * (M - 1) * (N - 1);
    unsigned int n2 = N * (M - 1) * (P - 1) + n1;
    unsigned int n3 = M * (N - 1) * (P - 1) + n2;
    typedef std::pair< unsigned int, short int > face_info;
    std::vector< face_info > all_faces(n3);
    unsigned int c = 0;
    for (unsigned int k = 0 ; k < P ; k++)
        for (unsigned int j = 0 ; j < N - 1 ; j++)
            for (unsigned int i = 0 ; i < M - 1 ; i++, c++) {
                unsigned int id = i + M * (j + N * k);
                all_faces[c] = face_info(id, 2);
            }
    for (unsigned int j = 0 ; j < N ; j++)
        for (unsigned int k = 0 ; k < P - 1 ; k++)
            for (unsigned int i = 0 ; i < M - 1 ; i++, c++) {
                unsigned int id = i + M * (j + N * k);
                all_faces[c] = face_info(id, 1);
            }
    for (unsigned int i = 0 ; i < M ; i++)
        for (unsigned int k = 0 ; k < P - 1 ; k++)
            for (unsigned int j = 0 ; j < N - 1 ; j++, c++) {
                unsigned int id = i + M * (j + N * k);
                all_faces[c] = face_info(id, 0);
            }

    int prev_pct = -1;
    std::vector< bool > conf_ok(M*N*P, false);
    std::vector< bool > ok_vertex(M*N*P, false);
    fixing_voxel = false;

    unsigned int nb_ok_vertices;
    if (read_info)
        nb_ok_vertices = read_vertex_info(ok_vertex, conf_ok);
    else
        nb_ok_vertices = compute_vertex_info(sample_grid, ok_vertex, conf_ok);
    std::cout << "there are " << nb_ok_vertices << " valid vertices ("
              << 100*nb_ok_vertices / (M*N*P) << "%)" << std::endl;
    std::cout << std::flush;

    unsigned int prev, cur;

    // -----------------------------------------------------------------------
    //
    // loop over all faces
    //
    // -----------------------------------------------------------------------

    crease::speedup = true;
    prev_pct = -1;

    // some timing material
    double total_empty_face_processing_time = 0;
    double total_found_face_processing_time = 0;
    unsigned int nb_found_face = 0;
    unsigned int nb_empty_face = 0;
    xavier::ProgressDisplay progress(false);
    progress.start();

    for (unsigned int f = 0 ; f < n3 ; f++) {
        FaceId fid = face_id(all_faces[f].first, all_faces[f].second, M, N);
        unsigned int i, j, k, axis;
        ijk(i, j, k, all_faces[f].first, M, N);
        axis = all_faces[f].second;
        std::string label;
        if (f < n1) label = "Z = cst";
        else if (f < n2) label = "Y = cst";
        else label = "X = cst";

        int pct = 100 * f / n3;
        if (pct > prev_pct) {
            std::cout << std::endl << "progress: " << pct << "% complete"
                      << std::endl << std::endl;
            prev_pct = pct;
        }

        // quick test first: is there at least one valid vertex?
        bool should_be_processed = false;
        for (unsigned int l = 0 ; l < 4 ; ++l) {
            unsigned int x, y, z;
            x = i + voxel_info::face_pt[axis][l][0];
            y = j + voxel_info::face_pt[axis][l][1];
            z = k + voxel_info::face_pt[axis][l][2];
            v[l] = sample_grid(x, y, z);
            // if (pct > 12 && !l) std::cerr << "v[0](" << x << ", " << y << ", " << z << ") = " << v[l] << '\n';
            ids[l] = x + M * (y + N * z);
            if (!conf_ok[ids[l]]) {
                // we skip faces that contain low confidence values
                should_be_processed = false;
                break;
            }
            else if (ok_vertex[ids[l]]) {
                should_be_processed = true;
            }
        }
        if (!should_be_processed) continue;

        for (unsigned int l = 0 ; l < 4 ; l++)
            ok_faces.push_back(v[l]);

        // do the actual work
        current_face_id = fid;
        std::vector< nvis::vec3 > points;
        current_face_info.set_info(v[0], v[1], v[2], v[3]);
        bool found_something;

        // std::cerr << "current valid face = " << v[0 ] << " - " << v[2] << '\n';

        progress.stop();
        progress.start();
        int res = search_face(points, v[0], v[1], v[2], v[3], max_depth, found_something);
        double dt = progress.cpu_time();

        if (res) {
            total_found_face_processing_time += dt;
            ++nb_found_face;
        }
        else {
            total_empty_face_processing_time += dt;
            ++nb_empty_face;
        }

        if (res) {
            std::cout << "extractor: PVO successful in " << dt << "sec. on face " << fid
                      << std::endl;
            nvis::vec3 refpos = points.front();
            double refval = the_wrapper->value(points.front());
            if (points.size() > 1) {
                ++nb_several;

                // we just select the position with the largest value.
                // an alternative possibility would be to select the position
                // associated with the largest crease strength.
                // note that either filtering strategy is discarding crease
                // lines that cross the same face twice
                for (unsigned int j = 1 ; j < points.size() ; j++) {
                    double val = crease::the_wrapper->value(points[j]);
                    if ((crease::is_ridge ? (val > refval) : (val < refval))) {
                        refval = val;
                        refpos = points[j];
                    }
                }
            }
            // add new entry to central database
            std::cout << "found a crease point: " << label << std::endl;
            all_face_points.push_back(refpos);
            face_point_ids[fid] = all_face_points.size() - 1;

            // add vertices involved in subdivision
            crease::vertices.push_back(current_vertices);

            crease::all_points_on_face.push_back(crease::point_on_face());
            crease::point_on_face& pof = crease::all_points_on_face.back();
            for (unsigned int n = 0 ; n < 4 ; n++) {
                pof.p[n] = v[n];
            }
            pof.e = the_wrapper->eigenvector(refpos, crease::is_ridge ? 0 : 2);
            pof.g = the_wrapper->gradient(refpos);
            pof.q = refpos;
        }
        else if (found_something) {
            face_point_ids[fid] = invalid_face_id; // don't look for it in future iterations
        }
    }
    std::cout << std::endl
              << "average time spent on face containing crease point is " << total_found_face_processing_time / (float)nb_found_face
              << std::endl
              << "average time spent on empty face is " << total_empty_face_processing_time / (float)nb_empty_face
              << std::endl;


    // -----------------------------------------------------------------------
    //
    // now we switch from a face-based to a voxel-based perspective to detect
    // existing loose ends
    //
    // -----------------------------------------------------------------------

    // loop over all voxels
    const unsigned int nb_voxels = (M - 1) * (N - 1) * (P - 1);
    apply_filter = false; // we are only going to filter the crease points not the vertices
    crease::speedup = false;

    std::vector< bool > is_fixed(nb_voxels, false);

    // store tangent vectors of crease points as induced by their cell neighbor
    std::map< FaceId, nvis::vec3 > tangents;

    int delta_x, delta_y, delta_z;
    delta_x = 1;
    delta_y = M - 1;
    delta_z = (M - 1) * (N - 1);
    const int f_neigh[6] = { -delta_z, delta_z, -delta_y, delta_x, delta_y, -delta_x };

    // following list keeps track of voxels containing loose ends
    std::list< unsigned int > voxels_to_look_into;

    // initially we check the status of all voxels and only process those
    // that contain two crease points, in which case we compute their tangent
    unsigned int nb_good = 0, nb_bad = 0;
    for (unsigned int v = 0 ; v < (M - 1)*(N - 1)*(P - 1) ; v++) {

        unsigned int k = v / (M - 1) / (N - 1);
        unsigned int j = (v - k * (M - 1) * (N - 1)) / (M - 1);
        unsigned int i = v % (M - 1);

        // check if voxel contains low confidence value
        if (!check_confidence(i, j, k, M, N, conf_ok)) {
            // std::cout << "skipping voxel with low confidence value" << std::endl;
            continue;
        }

        std::vector< unsigned int > fpids;
        std::vector< unsigned int > point_face;
        std::vector< FaceId > face_ids;

        bool found_invalid_solution = false;
        for (unsigned int f = 0 ; f < 6 ; f++) {
            FaceId fid = global_face_id(i, j, k, sample_grid, f);
            std::map< FaceId, unsigned int >::const_iterator it = face_point_ids.find(fid);
            if (it != face_point_ids.end()) {
                if (it->second != invalid_face_id) {
                    fpids.push_back(it->second);
                    point_face.push_back(f);
                    face_ids.push_back(fid);
                }
                else {
                    found_invalid_solution = true;
                }
            }
        }

        unsigned int np = fpids.size();
        if (np == 0) {
            is_fixed[v] = true;
            continue;
        }
        else if (np == 1 && !found_invalid_solution) {
            std::cout << "voxel #" << v << " will be processed later" << std::endl;
            voxels_to_look_into.push_back(v);
            ++nb_bad;
            continue;
        }
        else if (np == 2) {
            ++nb_good;
            std::cout << "segment found in voxel #" << v << std::endl;
            nvis::vec3 p0 = all_face_points[fpids[0]];
            nvis::vec3 p1 = all_face_points[fpids[1]];

            // debug code begins...
            round1.push_back(p0);
            round1.push_back(p1);
            // ...debug code ends

            nvis::vec3 dpdt = p0 - p1;
            dpdt /= nvis::norm(dpdt);
            add_segment(fpids);
            tangents[face_ids[0]] = dpdt;
            tangents[face_ids[1]] = -1 * dpdt; // orientation points to next voxel

            is_fixed[v] = true;
        }
        else if (np > 2) {
            ++nb_good;
            // connect points with highest (resp. lowest) value
            std::vector< double > vals(np);
            for (unsigned int n = 0 ; n < np ; ++n) {
                vals[n] = the_wrapper->value(all_face_points[fpids[n]]);
            }
            std::vector< unsigned int > sorted(np);
            xavier::sort(vals, sorted);
            unsigned int id0, id1;
            if (is_ridge) {
                id0 = sorted[np-1];
                id1 = sorted[np-2];
            }
            else {
                id0 = sorted[0];
                id1 = sorted[1];
            }

            std::cout << "TRICKY segment found in voxel #" << v << std::endl;
            nvis::vec3 p0 = all_face_points[fpids[id0]];
            nvis::vec3 p1 = all_face_points[fpids[id1]];

            // debug code begins...
            round1.push_back(p0);
            round1.push_back(p1);
            // ...debug code ends

            nvis::vec3 dpdt = p0 - p1;
            dpdt /= nvis::norm(dpdt);
            add_segment(fpids);
            tangents[face_ids[id0]] = dpdt;
            tangents[face_ids[id1]] = -1 * dpdt; // orientation points to next voxel

            is_fixed[v] = true;
        }
    }

    std::cout << "after initial inspection there are " << nb_good << " voxels containing a segment vs. "
              << nb_bad << " voxels containing a single point (" << 100.*(float)nb_good / (float)(nb_good + nb_bad)
              << "%)" << std::endl;


    // some timing material
    double total_failed_voxel_processing_time = 0;
    double total_fixed_voxel_processing_time = 0;
    double total_setup_time = 0;
    unsigned int nb_fixed_voxel = 0;
    unsigned int nb_failed_voxel = 0;
    unsigned int nb_setup = 0;

    // loop over all the voxels containing only one crease point
    {
        // keep track of which faces have undergone scrutiny.// this prevents us from doing the job twice if first attempt
        // was unsuccessful
        std::set< FaceId > scrutinized;

        unsigned int nb_proc = 0;
        unsigned int initial_nb = voxels_to_look_into.size();
        unsigned int nbp = 0;
        unsigned int list_sz = initial_nb;
        for (std::list< unsigned int >::iterator it = voxels_to_look_into.begin() ;
             it != voxels_to_look_into.end() ; it++, nbp++) {

            std::cout << "processing pending voxel #" << *it << " (" << nbp << ")"
                      << std::endl;

            the_wrapper->turn_buffer_on();

            ++nb_proc;

            unsigned int v = *it;
            unsigned int k = v / (M - 1) / (N - 1);
            unsigned int j = (v - k * (M - 1) * (N - 1)) / (M - 1);
            unsigned int i = v % (M - 1);

            std::vector< unsigned int > fpids;
            std::vector< unsigned int > point_face;
            std::vector< FaceId > face_ids;

            bool found_invalid_face = false;

            if (display_debug_info)
                std::cout << "looping over faces..." << std::endl;
            for (unsigned int f = 0 ; f < 6 ; f++) {
                FaceId fid = global_face_id(i, j, k, sample_grid, f);
                std::map< FaceId, unsigned int >::const_iterator it = face_point_ids.find(fid);
                if (it != face_point_ids.end()) {
                    if (it->second != invalid_face_id) {
                        if (display_debug_info)
                            std::cout << fid << " contains pos #" << it->second << std::endl;
                        fpids.push_back(it->second);
                        point_face.push_back(f);
                        face_ids.push_back(fid);
                    }
                    else {
                        found_invalid_face = true;
                    }
                }
                else {
                    if (display_debug_info)
                        std::cout << fid << " is empty" << std::endl;
                }
            }

            fixing_voxel = false;
            unsigned int np = fpids.size();
            if (np == 0) {
                // that cannot happen, by construction. let the user know that s/he has switched
                // to the twilight zone
                if (display_debug_info)
                    std::cout << "reached an empty voxel in second voxel loop!" << std::endl;
                is_fixed[v] = true;

                the_wrapper->turn_buffer_off();
                continue;
            }
            else if (np == 2) {
                // this case corresponds to a secondary fix/addition of a voxel as a side effect
                // of fixing its neighbor(s): we have killed two birds with one stone!
                std::cout << "segment found in voxel #" << v << std::endl;
                nvis::vec3 p0 = all_face_points[fpids[0]];
                nvis::vec3 p1 = all_face_points[fpids[1]];

                // debug code begins...
                round12.push_back(p0);
                round12.push_back(p1);
                // ...debug code ends

                nvis::vec3 dpdt = p0 - p1;
                dpdt /= nvis::norm(dpdt);
                add_segment(fpids);
                tangents[face_ids[0]] = dpdt;
                tangents[face_ids[1]] = -1. * dpdt; // orientation points to next voxel
                is_fixed[v] = true;
            }
            else if (np > 2) {
                // connect points with highest (resp. lowest) value
                std::vector< double > vals(np);
                for (unsigned int n = 0 ; n < np ; ++n) {
                    vals[n] = the_wrapper->value(all_face_points[fpids[n]]);
                }
                std::vector< unsigned int > sorted(np);
                xavier::sort(vals, sorted);
                unsigned int id0, id1;
                if (is_ridge) {
                    id0 = sorted[np-1];
                    id1 = sorted[np-2];
                }
                else {
                    id0 = sorted[0];
                    id1 = sorted[1];
                }

                std::cout << "TRICKY segment found in voxel #" << v << std::endl;
                nvis::vec3 p0 = all_face_points[fpids[id0]];
                nvis::vec3 p1 = all_face_points[fpids[id1]];

                // debug code begins...
                round1.push_back(p0);
                round1.push_back(p1);
                // ...debug code ends

                nvis::vec3 dpdt = p0 - p1;
                dpdt /= nvis::norm(dpdt);
                add_segment(fpids);
                tangents[face_ids[id0]] = dpdt;
                tangents[face_ids[id1]] = -1 * dpdt; // orientation points to next voxel

                is_fixed[v] = true;
            }
            else if (np == 1) {
                if (found_invalid_face) {
                    is_fixed[v] = true;
                    continue;
                }

                // debug code begins...
                if (display_debug_info)
                    std::cout << std::endl << "must now fix voxel #" << v << std::endl;
                // ...debug code ends

                bool fixed = false;

                progress.stop();
                progress.start();

                // build up voxel geometry
                std::vector< nvis::vec3 > vert(8);
                sample_grid.voxel(vert, i, j, k);
                if (display_debug_info)
                    std::cout << "voxel is: 0=" << vert[0] << ", 1=" << vert[1] << ", 2=" << vert[2]
                              << ", 3=" << vert[3] << ", 4=" << vert[4] << ", 5=" << vert[5] << ", 6=" << vert[6]
                              << ", 6=" << vert[6] << ", 7=" << vert[7] << std::endl;

                fixing_voxel = true;

                std::pair< int, nvis::vec3 > first_guess; // face id / position on face
                first_guess.first = -1;   // no first guess so far

                // check if we have a tangent available at this position
                std::map< FaceId, nvis::vec3 >::iterator tang_it = tangents.find(face_ids[0]);

                // tracking solution specific to crease lines of mode
                double h = 0.25;
                if (crease_kind == 2) {
                    unsigned int fid_out;
                    nvis::vec3 out;
                    if (track_ridge_line(out, fid_out, vert, all_face_points[fpids[0]], point_face[0])) {
                        double val = the_wrapper->value(out);
                        if (display_debug_info)
                            std::cout << "tracking was successful. end position = " << out
                                      << ", face id = " << fid_out << ", value = " << val << std::endl;
                        first_guess.first = fid_out;
                        first_guess.second = out;

                        // tracking approach provides us with greater accuracy
                        h = 0.1;
                    }
                    else {
                        if (display_debug_info)
                            std::cout << "tracking failed" << std::endl;
                    }
                }
                else if (tang_it != tangents.end()) {
                    unsigned int fid_out;
                    nvis::vec3 out;
                    trace_ray(out, fid_out, vert, all_face_points[fpids[0]], tang_it->second, point_face[0]);
                    double val = the_wrapper->value(out);
                    if (display_debug_info)
                        std::cout << "ray tracing was successful. end position = " << out
                                  << ", face id = " << fid_out << ", value = " << val << std::endl;
                    first_guess.first = fid_out;
                    first_guess.second = out;
                }

                bool found = false;
                bool found_something = false;
                std::vector< nvis::vec3 > points;
                pvo_solution the_solution;

                // if we choose to trust the tracking outcome we take it as solution
                if (crease_kind == 2 && bold_move && first_guess.first >= 0) {
                    double val = the_wrapper->value(first_guess.second);
                    if ((is_ridge ? (val > value_threshold_select) :
                         (val < value_threshold_select))) {
                        found = true;
                        the_solution.p = first_guess.second;
                        the_solution.local_fid = first_guess.first;
                        the_solution.global_fid = global_face_id(i, j, k, sample_grid, first_guess.first);
                    }
                }

                progress.stop();
                total_setup_time += progress.cpu_time();
                ++nb_setup;

                progress.start();
                // if we have identified a first guess, look for it in priority
                if (!found && first_guess.first >= 0) {
                    fixing_voxel = true;
                    if (display_debug_info)
                        std::cout << "we are about to scrutinize the right face" << std::endl;
                    const unsigned int* ids = faces[first_guess.first];
                    face_type face(vert[ids[0]], vert[ids[1]], vert[ids[2]], vert[ids[3]]);

                    if (crease_kind == 2) {
                        face_type out;
                        refine_face(out, face, first_guess.second, h);
                        if (display_debug_info)
                            std::cout << std::endl << std::endl << "we are scrutinizing face: "
                                      << out.p[0] << ", " << out.p[1] << ", " << out.p[2] << ", " << out.p[3]
                                      << std::endl;

                        unsigned int n = 0;
                        if (crease_kind < 2) n = 2;
                        for (; n < 3 && !found ; ++n) {

                            if (!n) {
                                if (display_debug_info)
                                    std::cout << "first attempt: bold approach" << std::endl;
                                xavier::crease::speedup = true;
                            }
                            else if (n == 1) {
                                if (display_debug_info)
                                    std::cout << "first attempt failed: cautious approach on second attempt" << std::endl;
                                xavier::crease::speedup = false;
                            }
                            else {
                                if (display_debug_info)
                                    std::cout << "everything we tried so far failed. using a smaller area" << std::endl;
                                xavier::crease::speedup = false;
                                h *= 0.5;
                                refine_face(out, face, first_guess.second, h);
                            }
                            current_vertices.clear();
                            current_face_info.set_info(out.p[0], out.p[1], out.p[2], out.p[3]);
                            found = search_face(points, out.p[0], out.p[1], out.p[2], out.p[3], max_depth, found_something);

                            if (found) {
                                the_solution.p = points[0];
                                the_solution.local_fid = first_guess.first;
                                the_solution.global_fid = global_face_id(i, j, k, sample_grid, first_guess.first);
                            }
                        }
                    }
                    else { // simply look on face pointed to by tangent vector
                        if (display_debug_info)
                            std::cout << std::endl << std::endl << "we are scrutinizing face: "
                                      << face.p[0] << ", " << face.p[1] << ", " << face.p[2] << ", " << face.p[3]
                                      << std::endl;

                        xavier::crease::speedup = false;
                        current_vertices.clear();
                        current_face_info.set_info(face.p[0], face.p[1], face.p[2], face.p[3]);
                        found = search_face(points, face.p[0], face.p[1], face.p[2], face.p[3], max_depth_fix,
                                            found_something);

                        FaceId _fid = global_face_id(i, j, k, sample_grid, first_guess.first);
                        if (found) {
                            the_solution.p = points[0];
                            the_solution.local_fid = first_guess.first;
                            the_solution.global_fid = _fid;
                        }
                        else {
                            scrutinized.insert(_fid);
                        }
                    }
                }

                // brute force loop over all faces, if necessary
                unsigned int __d;
                // for (__d = 1 ; __d <= max_depth_fix && !found && !found_something ; ++__d, crease::max_int_error /= 2) {
                {
                    __d = max_depth_fix;
                    crease::max_int_error /= 5;
                    xavier::crease::speedup = false;
                    xavier::crease::fixing_voxel = true;
                    for (unsigned int vf = 0 ; vf < 6 && !found && !found_something; vf++) {
                        if (vf == point_face[0] || // skip position that we already know
                            (crease_kind < 2 && // or face we have already ruled out
                             first_guess.first >= 0 && vf == first_guess.first))
                            continue;
                        fixing_voxel = true;
                        FaceId fid = global_face_id(i, j, k, sample_grid, vf);
                        if (scrutinized.find(fid) != scrutinized.end())
                            continue;
                        const unsigned int* face = voxel_info::faces[vf];

                        if (display_debug_info)
                            std::cout << std::endl << std::endl << "we are scrutinizing face: "
                                      << vert[face[0]] << ", " << vert[face[1]] << ", "
                                      << vert[face[2]] << ", " << vert[face[3]] << std::endl << std::endl;

                        current_vertices.clear();
                        current_face_info.set_info(vert[face[0]], vert[face[1]], vert[face[2]], vert[face[3]]);
                        found = search_face(points, vert[face[0]], vert[face[1]],
                                            vert[face[2]], vert[face[3]], __d, found_something);
                        if (found) {
                            the_solution.p = points[0];
                            the_solution.local_fid = vf;
                            the_solution.global_fid = fid;
                        }
                        else if (found_something || __d == max_depth_fix) {
                            // we have reached the maximum allowed subdivision depth
                            // without being successful. we won't be any more successful
                            // in further attempts
                            scrutinized.insert(fid);
                        }
                    }
                }

                the_wrapper->turn_buffer_off();

                progress.stop();
                double dt = progress.cpu_time();

                if (found) {
                    // add new entry to central database
                    std::cout << "found a crease point that was hiding #"
                              << added_vertices.size() + 1 << "! : " << the_solution.p << std::endl;
                    std::cout << "it took a subdivision depth = " << __d << std::endl;
                    all_face_points.push_back(the_solution.p);
                    face_point_ids[the_solution.global_fid] = all_face_points.size() - 1;
                    std::cout << "we added pos #" << all_face_points.size() - 1 << " to face " << the_solution.global_fid << std::endl;
                    added_vertices.push_back(all_face_points.size() - 1);
                    fpids.push_back(all_face_points.size() - 1);
                    point_face.push_back(the_solution.local_fid);
                    face_ids.push_back(the_solution.global_fid);

                    // update tangent information
                    nvis::vec3 p0 = all_face_points[fpids[0]];
                    nvis::vec3 p1 = all_face_points[fpids[1]];

                    // debug code begins...
                    round2.push_back(p0);
                    round2.push_back(p1);
                    // ...debug code ends

                    nvis::vec3 dpdt = p0 - p1;
                    dpdt /= nvis::norm(dpdt);
                    add_segment(fpids);
                    tangents[face_ids[0]] = dpdt;
                    tangents[face_ids[1]] = -1 * dpdt; // orientation points to next voxel
                    is_fixed[v] = true;

                    int modif = v + f_neigh[the_solution.local_fid];
                    if (is_fixed[modif]) {
                        voxels_to_look_into.push_back(modif);
                        ++list_sz;
                        // debug code begins...
                        if (display_debug_info)
                            std::cout << "we will have to look further into voxel #"
                                      << modif << " (" << nb_voxels << ")" << std::endl;
                        // ...debug code ends
                    }
                    std::cout << nb_proc << " voxels processed so far, " << (list_sz - nb_proc) << " voxels to go" << std::endl;

                    // add vertices involved in subdivision
                    crease::vertices.push_back(current_vertices);

                    // add voxel to list of fixed ones for display
                    for (unsigned int l = 0 ; l < 8 ; l++) {
                        fixed_voxels.push_back(vert[l]);
                    }

                    fixed = true;
                }
                else if (found_something) {

                    // nothing to be done for this voxel where a solution exists, as required by
                    // the continuity of the crease lines, but does not meet our filtering criteria
                }
                else  {
                    if (display_debug_info) {
                        std::cout << "we failed to fix this voxel" << std::endl;
                        std::cout << "vertices are: " << std::endl;
                    }
                    std::cout << "FAILED" << std::endl;
                    for (unsigned int l = 0 ; l < 8 ; l++) {
                        problematic_voxels.push_back(vert[l]);
                        if (display_debug_info)
                            std::cout << "#" << l << ": " << vert[l] << std::endl;
                    }

                    // failed_voxels.push_back(v);
                }

                if (found || found_something) {
                    total_fixed_voxel_processing_time += dt;
                    ++nb_fixed_voxel;
                }
                else {
                    total_failed_voxel_processing_time += dt;
                    ++nb_failed_voxel;
                }
            }
        }
        // voxels_to_look_into.swap(failed_voxels);
        //     failed_voxels.clear();
    }
    std::cout << "average time spent on fixed voxels = " << total_fixed_voxel_processing_time / (float)nb_fixed_voxel
              << " sec." << std::endl
              << "average time spent on failed voxels = " << total_failed_voxel_processing_time / (float)nb_failed_voxel
              << " sec." << std::endl;

    // identify connected components
    connect_segments(creases);

    unsigned int nb_cells = (M - 1) * (N - 1) * (P - 1);
    std::cout << "percentage of candidate voxels: "
              << 100*(double)nbok / (double)nb_cells
              << std::endl;
    std::cout << "percentage of voxels discarded because of value: "
              << 100*(double)nb_wrong_val / (double)n3 << std::endl
              << "percentage of voxels discarded because of strength: "
              << 100*(double)nb_wrong_eval / (double)n3 << std::endl
              << "number of failed convergences: " << failed_conv << std::endl
              << "percentage of faces containing several zero crossing: "
              << 100*(double)nb_several / (double)all_face_points.size() << std::endl
              << "number of PVO computations performed: " << xavier::crease::nb_pvo << std::endl;

    std::cout << "number of segments = " << all_edges.size()
              << ", number of connected components: " << creases.size()
              << std::endl;
}

// --------------------------------------------------------------------------

void crease::connect_segments(std::vector< line >& creases)
{
// compute point <- point -> point connections
    std::vector< std::pair< int, int > > connected_to(all_face_points.size(),
            std::pair< int, int >(-1, -1));
    for (unsigned int i = 0 ; i < all_edges.size() ; i++) {
// face indices for current edge
        unsigned int f0, f1;
        f0 = all_edges[i].first;
        f1 = all_edges[i].second;

        int a, b;
        a = connected_to[f0].first;
        b = connected_to[f1].first;

        if (a == -1)
            connected_to[f0].first = f1;
        else
            connected_to[f0].second = f1;

        if (b == -1)
            connected_to[f1].first = f0;
        else
            connected_to[f1].second = f0;
    }

    creases.clear();
    std::vector< bool > inserted(all_face_points.size(), false);
    for (unsigned int i = 0 ; i < all_edges.size() ; i++) {
// edge end points
        unsigned int i0, i1;
        i0 = all_edges[i].first;
        i1 = all_edges[i].second;

        assert(i0 < all_face_points.size() && i1 < all_face_points.size());

        if (inserted[i0] || inserted[i1]) continue;

        unsigned int cur, prev;
        int link0, link1;

// start a new connected component
        std::list< unsigned int > my_list;

// initialize connected component with these two points
        my_list.push_back(i0);
        my_list.push_back(i1);

// append forward
        prev = i0;
        cur = i1;

        inserted[prev] = true;
        for (; ;) {
            inserted[cur] = true;
            link0 = connected_to[cur].first; // always >= 0
            link1 = connected_to[cur].second;
            assert(link0 >= 0);
            if ((unsigned int)link0 != prev && !inserted[link0])
                my_list.push_back(link0);
            else if (link1 >= 0 && (unsigned int)link1 != prev && !inserted[link1])
                my_list.push_back(link1);
            else break;

            prev = cur;
            cur = my_list.back();
        }

// append backward
        cur = i0;
        prev = i1;
        for (; ;) {
            inserted[cur] = true;
            link0 = connected_to[cur].first;
            link1 = connected_to[cur].second;
            assert(link1 < 0 || link1 < all_face_points.size());
            if ((unsigned int)link0 != prev && !inserted[link0])
                my_list.push_front(link0);
            else if (link1 >= 0 && (unsigned int)link1 != prev && !inserted[link1])
                my_list.push_front(link1);
            else break;

            prev = cur;
            cur = my_list.front();
        }

        creases.push_back(std::vector< unsigned int >(my_list.begin(), my_list.end()));
    }

    std::cout << "verifying results:" << std::endl;
    for (unsigned int i = 0 ; i < creases.size() ; i++) {
        std::cout << "component #" << i << ": (" << creases[i].size() << ")" << std::flush;
        for (unsigned int j = 0 ; j < creases[i].size() ; j++) {
            std::cout << creases[i][j] << " " << std::flush;
        }
        std::cout << std::endl;
    }
}

// --------------------------------------------------------------------------
