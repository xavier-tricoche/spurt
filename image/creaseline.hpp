#ifndef __CREASE_LINE_HPP__
#define __CREASE_LINE_HPP__

#include <ext/spurt/image/probe.hpp>
#include <ext/nvis/math/fixed_vector.hpp>
#include <teem/nrrd.h>
#include <set>
#include <map>
#include <vector>

namespace spurt {
namespace crease {
std::vector< std::pair< nvis::vec3, nvis::vec3 > > problematic_edges;
std::vector< nvis::vec3 > face_points;
std::vector< std::pair< unsigned int, unsigned int > > edges;
double threshold;
unsigned int subdiv;

struct FaceId {
    FaceId(unsigned int i0, unsigned int i1,
           unsigned int i2, unsigned int i3)
            : is(4) {
        is[0] = i0;
        is[1] = i1;
        is[2] = i2;
        is[3] = i3;
        std::sort(is.begin(), is.end());
    }

    FaceId(const FaceId& fid)
            : is(4) {
        for (unsigned int i = 0 ; i < 4 ; i++)
            is[i] = fid.is[i];
    }

    int operator<(const FaceId& fid) const {
        return (is[0] < fid.is[0] ||
                (is[0] == fid.is[0] &&
                 (is[1] < fid.is[1] ||
                  (is[1] == fid.is[1] &&
                   (is[2] < fid.is[2] ||
                    (is[2] == fid.is[2] &&
                     is[3] < fid.is[3]))))));
    }

    std::vector< unsigned int > is;
};

ostream& operator<< (ostream& os, const FaceId& fid)
{
    os << "[ " << fid.is[0] << ", " << fid.is[1] << ", "
    << fid.is[2] << ", " << fid.is[3] << "]";
    return os;
}

std::vector< nvis::vec3 > non_orientable;

struct Grid {
    Grid(unsigned int size[3], double min[3], double max[3]) {
        for (unsigned int i = 0 ; i < 3 ; i++) {
            if (false) {
                // physical coordinates
                _min[i] = min[i];
                _max[i] = max[i];
                _size[i] = size[i];
                _d[i] = (_max[i] - _min[i]) / (double)(_size[i] - 1);
            }
            else {
                // raster coordinates
                _min[i] = 0;
                _max[i] = size[i] - 1;
                _size[i] = size[i];
                _d[i] = 1;
            }
        }
    }

    unsigned int id(unsigned int i, unsigned int j, unsigned int k) {
        return i + _size[0]*(j + _size[1]*k);
    }

    nvis::vec3 operator()(unsigned int i, unsigned int j, unsigned int k) {
        return nvis::vec3(_min[0] + i*_d[0],
                          _min[1] + j*_d[1],
                          _min[2] + k*_d[2]);
    }

    double _min[3], _max[3];
    unsigned int _size[3];
    double _d[3];
};

bool zero_crossing(std::vector< nvis::vec3 >& zeros,
                   const nvis::vec3& p0, const nvis::vec3& p1,
                   const spurt::gage_interface::scalar_wrapper& gH,
                   unsigned int idx, unsigned int alt_idx);
                   
void zero_crossing_e13(std::vector < std::pair < nvis::vec3,
                       nvis::vec3 > > & zeros,
                       const nvis::vec3& p0, const nvis::vec3& p1,
                       const spurt::gage_interface::scalar_wrapper& gH,
                       unsigned int idx, unsigned int ref);
                       
void zero_crossing_e2(std::vector< nvis::vec3 >& zeros,
                      const nvis::vec3& p0, const nvis::vec3& p1,
                      const nvis::vec3& ev0, const nvis::vec3& ev1,
                      const spurt::gage_interface::scalar_wrapper& gH,
                      unsigned int ref);
                      
void edge_orientation(nvis::vec3& ev1, const nvis::vec3& reference,
                      const nvis::vec3& p0, const nvis::vec3& p1,
                      const spurt::gage_interface::scalar_wrapper& gH,
                      unsigned int idx);
                      
bool face_orientation(std::vector< nvis::vec3 >& vecs,
                      const std::vector< nvis::vec3 >& p,
                      const spurt::gage_interface::scalar_wrapper& gH,
                      unsigned int idx, bool ridge);
                      
void extract_lines(const Nrrd* nrrd, bool ridge);

};
};

void spurt::crease::
edge_orientation(nvis::vec3& ev1, const nvis::vec3& reference,
                 const nvis::vec3& p0, const nvis::vec3& p1,
                 const spurt::gage_interface::scalar_wrapper& gH,
                 unsigned int idx)
{
    using namespace spurt;
    using namespace gage_interface;
    using namespace nvis;
    using namespace std;

    vector< vec3 > evecs(3);
    vec3 p, prev, cur;
    prev = cur = reference;
    if (subdiv)
        for (unsigned int i = 1 ; i < subdiv + 2 ; i++) {
            double u = (double)i / (double)(subdiv + 1);
            p = (1 - u) * p0 + u * p1;
            gH.hess_evecs(p, evecs);
            cur = evecs[idx];
            if (inner(cur, prev) < 0)
                cur *= -1;
            prev = cur;
        }
    else {
        gH.hess_evecs(p1, evecs);
        cur = evecs[idx];
        if (inner(cur, prev) < 0)
            cur *= -1;
    }
    ev1 = cur;
}

bool spurt::crease::
face_orientation(std::vector< nvis::vec3 >& vecs,
                 const std::vector< nvis::vec3 >& p,
                 const spurt::gage_interface::scalar_wrapper& gH,
                 unsigned int idx,
                 bool ridge)
{
    using namespace spurt;
    using namespace gage_interface;
    using namespace nvis;
    using namespace std;

    vecs.resize(4);
    vector< vec3 > evecs(3);
    vec3 evals;
    gH.hess_evecs(p[0], evecs);
    gH.hess_evals(p[0], evals);
    if ((ridge && fabs(evals[1] - evals[2]) / fabs(evals[2]) < 0.01) ||
        (!ridge && fabs(evals[1] - evals[0]) / fabs(evals[0]) < 0.01))
        cout << "ambiguous eigenvalues!!" << endl;

    vecs[0] = evecs[idx];

// propagate orientation at 1st vertex to the others
    edge_orientation(vecs[1], vecs[0], p[0], p[1], gH, idx);
    edge_orientation(vecs[3], vecs[0], p[0], p[3], gH, idx);
    edge_orientation(vecs[2], vecs[1], p[1], p[2], gH, idx);

// check if found orientation is overall consistent
    vec3 test;
    edge_orientation(test, vecs[3], p[3], p[2], gH, idx);
    if (inner(test, vecs[2]) < 0) return false;

    return true;
}

void spurt::crease::
zero_crossing_e13(std::vector< std::pair< nvis::vec3, nvis::vec3 > >& zeros,
                  const nvis::vec3& p0, const nvis::vec3& p1,
                  const spurt::gage_interface::scalar_wrapper& gH,
                  unsigned int idx, unsigned int ref)
{
    using namespace spurt;
    using namespace gage_interface;
    using namespace nvis;
    using namespace std;

// create arbitrary basis for 2D eigenspace
    vector< vec3 > evecs(3);

// determine reference values at edge corners
    vec3 ev0, ev1;
    gH.hess_evecs(p0, evecs);
    ev0 = evecs[idx];
    gH.hess_evecs(p1, evecs);
    ev1 = evecs[idx];
    if (inner(ev0, ev1) < 0) {
        ev1 *= -1;
    }

    double t = 0;
    double dt = 1 / (double)subdiv;
    double dot0, dot1;
    vec3 grad;
    gH.gradient(p0, grad);
    dot0 = inner(grad, ev0);
    while (t < 1) {
        double tt = t + dt;
        if (tt > 1) tt = 1;
        vec3 p = (1 - tt) * p0 + tt * p1;
        gH.hess_evecs(p, evecs);
        // interpolate corner values
        vec3 ev = (1 - tt) * ev0 + tt * ev1;
        // project onto 2D eigenspace
        vec3 reference = evecs[ref];
        ev -= (ev * reference) * reference;
        ev /= norm(ev);

        gH.gradient(p, grad);
        dot1 = inner(grad, ev);
        if (dot0*dot1 < 0) {
            double u = -dot0 / (dot1 - dot0);
            double v = (1 - u) * t + u * tt;
            vec3 q = (1 - v) * p0 + v * p1;
            // compute eigenvector at corresponding location
            gH.hess_evecs(q, evecs);
            ev = (1 - v) * ev0 + v * ev1;
            ev -= (ev * evecs[ref]) * evecs[ref];
            ev /= norm(ev);
            zeros.push_back(pair< vec3, vec3 >(q, ev));
        }

        dot0 = dot1;
        t = tt;
    }
}

void spurt::crease::
zero_crossing_e2(std::vector< nvis::vec3 >& zeros,
                 const nvis::vec3& p0, const nvis::vec3& p1,
                 const nvis::vec3& ev0, const nvis::vec3& ev1,
                 const spurt::gage_interface::scalar_wrapper& gH,
                 unsigned int ref)
{
    using namespace spurt;
    using namespace gage_interface;
    using namespace nvis;
    using namespace std;

    vector< vec3 > evecs(3);
    vec3 grad, ev;
    gH.hess_evecs(p0, evecs);
    ev = cross(evecs[ref], ev0);

    double t = 0;
    double dt = 1 / (double)subdiv;

    gH.gradient(p0, grad);
    double dot0, dot1;
    dot0 = inner(grad, ev);

    while (t < 1) {
        double tt = t + dt;
        if (tt > 1) tt = 1;
        vec3 p = (1 - tt) * p0 + tt * p1;
        gH.hess_evecs(p, evecs);
        ev = (1 - tt) * ev0 + tt * ev1;
        ev -= (ev * evecs[ref]) * evecs[ref];
        ev /= norm(ev);   // fake "major" or "minor" eigenvector
        // obtain corresponding medium eigenvector by rotation in 2D eigenspace
        ev = cross(evecs[ref], ev);
        gH.gradient(p, grad);
        dot1 = inner(grad, ev);
        if (dot0*dot1 < 0) {
            double u = -dot0 / (dot1 - dot0);
            double v = (1 - u) * t + u * tt;
            vec3 q = (1 - v) * p0 + v * p1;
            zeros.push_back(q);
        }

        t = tt;
        dot0 = dot1;
    }
}


bool spurt::crease::
zero_crossing(std::vector< nvis::vec3 >& zeros,
              const nvis::vec3& p0, const nvis::vec3& p1,
              const spurt::gage_interface::scalar_wrapper& gH,
              unsigned int idx, unsigned int alt_idx)
{
    using namespace spurt;
    using namespace gage_interface;
    using namespace nvis;
    using namespace std;

    zeros.clear();

    cout << "looking for zero crossings along an edge" << endl;

    vector< vec3 > eigenvectors(3);
    vec3 gradient;
    vec3 prev, cur;
    gH.hess_evecs(p0, eigenvectors);
    prev = eigenvectors[idx];
    double dot0, dot1;
    gH.gradient(p0, gradient);
    dot0 = inner(gradient, prev);

    double t = 0;
    double dt = 0.5;
    while (t < 1) {
        double tt = t + dt;
        if (tt > 1) tt = 1;

        cout << "t = " << tt << endl;

        vec3 p = (1 - tt) * p0 + tt * p1;
        gH.hess_evecs(p, eigenvectors);
        cur = eigenvectors[idx];
        if (inner(cur, prev) < 0) {   // consistent orientation
            cur *= -1;
        }
        if (inner(cur, prev) < 0.85) {   // small angular variation
            dt *= 0.5;
            if (dt < 0.001) {
                cout << "step size underflow in zero_crossing" << endl;

                cout << "switching to a different eigenvector" << endl;
                cur = eigenvectors[alt_idx];
                if (inner(cur, prev) < 0) {
                    cur *= -1;
                }
                if (inner(cur, prev) < 0.85) {
                    cout << "switching to other eigenvector didn't help"
                    << endl;

                    problematic_edges.push_back(pair< vec3, vec3 >(p0, p1));
                    return false;
                }
                else {
                    // switch eigenvectors
                    unsigned int tmp = alt_idx;
                    alt_idx = idx;
                    idx = tmp;
                    dt = 0.5;
                }
            }
            continue;
        }

        gH.gradient(p, gradient);
        dot1 = inner(gradient, cur);
        if (dot1*dot0 < 0) {
            double u = -dot0 / (dot1 - dot0);
            double v = (1 - u) * t + u * tt;
            vec3 q = (1 - v) * p0 + v * p1;

            zeros.push_back(q);
            dt = 0.5;
        }

        t = tt;
        prev = cur;
        dot0 = dot1;
    }

    return true;
}

void spurt::crease::extract_lines(const Nrrd* nrrd, bool ridge)
{
    using namespace spurt;
    using namespace gage_interface;
    using namespace nvis;
    using namespace std;

// collect information about size and bounding box
    unsigned int size[3];
    double min[3], max[3];
    for (unsigned int i = 0 ; i < 3 ; i++) {
        min[i] = nrrd->axis[i].min;
        max[i] = nrrd->axis[i].max;
        size[i] = nrrd->axis[i].size;
    }
    Grid grid(size, min, max);
    scalar_wrapper gH(nrrd);

    set< FaceId > face_done;
    map< FaceId, unsigned int > face_found;
    face_points.clear();
    non_orientable.clear();

// loop over mesh cells
    unsigned int faces[6][4] = { { 0, 1, 2, 3 },
        { 4, 5, 6, 7 },
        { 0, 1, 5, 4 },
        { 1, 2, 6, 5 },
        { 2, 3, 7, 6 },
        { 3, 0, 4, 7 }
    };

    unsigned int nbok = 0;

// compute some statistics
    vec3 evals;
    std::vector< double > means;
//unsigned int counter=0;
    for (unsigned int i = 0 ; i < size[0] ; i++)
        for (unsigned int j = 0 ; j < size[1] ; j++)
            for (unsigned int k = 0 ; k < size[2] ; k++) {
                gH.hess_evals(grid(i, j, k), evals);
                if ((ridge && evals[1] < 0) ||
                    (!ridge && evals[1] > 0)) {
                    means.push_back(evals[1]);
                }
            }

    sort(means.begin(), means.end());
    unsigned int _id = (ridge
                        ? (unsigned int)(means.size() * threshold)
                        : (unsigned int)(means.size() * (1 - threshold)));
    double _threshold = means[_id];

    cout << "threshold = " << _threshold << endl;
    cout << "median valid value =  " << means[means.size()/2] << endl;

//unsigned int nb_valid = 0;
//unsigned int nb_wrong_evals = 0;
//unsigned int nb_wrong_espace = 0;
//unsigned int nb_wrong_basis = 0;

    for (unsigned int k = 0 ; k < size[2] - 1 ; k++)
        for (unsigned int j = 0 ; j < size[1] - 1 ; j++)
            for (unsigned int i = 0 ; i < size[0] - 1 ; i++) {
                // voxel vertices
                vector< vec3 > v(8);
                v[0] = grid(i  , j  , k);
                v[1] = grid(i + 1, j  , k);
                v[2] = grid(i + 1, j + 1, k);
                v[3] = grid(i  , j + 1, k);
                v[4] = grid(i  , j  , k + 1);
                v[5] = grid(i + 1, j  , k + 1);
                v[6] = grid(i + 1, j + 1, k + 1);
                v[7] = grid(i  , j + 1, k + 1);

                vector< unsigned int > ids(8);
                ids[0] = grid.id(i  , j  , k);
                ids[1] = grid.id(i + 1, j  , k);
                ids[2] = grid.id(i + 1, j + 1, k);
                ids[3] = grid.id(i  , j + 1, k);
                ids[4] = grid.id(i  , j  , k + 1);
                ids[5] = grid.id(i + 1, j  , k + 1);
                ids[6] = grid.id(i + 1, j + 1, k + 1);
                ids[7] = grid.id(i  , j + 1, k + 1);

                // check sign of eigenvalues
                vec3 evals;
                bool candidate = true;
                for (unsigned int l = 0 ; l < 8 ; l++) {
                    gH.hess_evals(v[l], evals);
                    if ((ridge && evals[1] >= _threshold) ||
                        (!ridge && evals[1] <= _threshold)) {
                        candidate = false;
                        break;
                    }
                }
                if (!candidate) continue;

//       std::cout << "valid voxel" << std::endl;

                nbok++;

                // loop over voxel faces
                vector< unsigned int > face_point;
                for (unsigned int f = 0 ; f < 6 ; f++) {
                    // face vertices
                    vector< vec3 > p(4);
                    for (unsigned int i = 0 ; i < 4 ; i++)
                        p[i] = v[faces[f][i]];

                    // check if we have processed that face already
                    FaceId fid(ids[faces[f][0]],
                               ids[faces[f][1]],
                               ids[faces[f][2]],
                               ids[faces[f][3]]);
                    if (face_done.find(fid) != face_done.end()) {
                        // add corresponding position, if any
                        map< FaceId, unsigned int >::const_iterator it =
                            face_found.find(fid);
                        if (it != face_found.end()) {
                            face_point.push_back(it->second);
                        }

                        // move to next voxel face
                        continue;
                    }

                    // mark current face as processed
                    face_done.insert(fid);

                    vector< pair< vec3, vec3 > > edge_point;

                    for (unsigned int e = 0 ; e < 4 ; e++) {
                        unsigned int i0 = e;
                        unsigned int i1 = (e + 1) % 4;
                        vector< pair< vec3, vec3 > > zeros_e13;
                        zero_crossing_e13(zeros_e13, p[i0], p[i1], gH,
                                          (ridge ? 2 : 0),
                                          (ridge ? 0 : 2));
                        if (zeros_e13.size() > 1) {
                            cout << "found several (" << zeros_e13.size()
                            << ") zero crossings along a single edge"
                            << endl;
                            continue;
                        }
                        if (zeros_e13.size() == 1)
                            edge_point.push_back(zeros_e13[0]);
                    }
                    if (!edge_point.size()) continue;
                    if (edge_point.size() != 2) {
                        cout << "found " << edge_point.size()
                        << " point(s) on face edges"
                        << endl;
                        continue;
                    }

                    // orient medium eigenvector along found edge

                    /*
                        vector< vec3 > vecs(3);
                        gH.hess_evecs( edge_point[0], vecs );
                        vec3 ev0 = vecs[1];
                        vec3 ev1;
                        edge_orientation( ev1, ev0, edge_point[0], edge_point[1],
                            gH, 1 );

                    // check for zero crossing of <grad,e>
                        vec3 g0, g1;
                        gH.gradient( edge_point[0], g0 );
                        gH.gradient( edge_point[1], g1 );
                        double v0 = inner( g0, ev0 );
                        double v1 = inner( g1, ev1 );
                        if ( v0*v1 < 0 )
                        {
                            double u = -v0/(v1-v0);
                            vec3 q = (1-u)*edge_point[0] + u*edge_point[1];
                            face_points.push_back( q );
                            face_point.push_back( face_points.size()-1 );
                            face_found[fid] = face_points.size()-1;
                        }
                    */

                    vector< vec3 > zeros;
                    zero_crossing_e2(zeros,
                                     edge_point[0].first,
                                     edge_point[1].first,
                                     edge_point[0].second,
                                     edge_point[1].second,
                                     gH, (ridge ? 0 : 2));
                    if (zeros.size() > 1) {
                        cout << "found " << zeros.size() << " crease points on "
                        << "voxel face" << endl;
                        continue;
                    }
                    if (zeros.size() == 1) {
                        face_points.push_back(zeros[0]);
                        face_point.push_back(face_points.size() - 1);
                        face_found[fid] = face_points.size() - 1;
                    }
                }

                if (face_point.size() &&
                    face_point.size() != 2) {
                    cout << "found " << face_point.size()
                    << " crease point(s) on voxel faces"
                    << endl;
                }
                else if (face_point.size() == 2) {
                    cout << "found a crease line!!" << endl;
                    edges.push_back(pair< unsigned int, unsigned int >
                                    (face_point[0], face_point[1]));
                    /*
                        output << "p "
                            << face_point[0][0] << " "
                            << face_point[0][1] << " "
                            << face_point[0][2] << endl
                            << "p "
                            << face_point[1][0] << " "
                            << face_point[1][1] << " "
                            << face_point[1][2] << endl
                            << "n" << endl;
                    */

                }
            }

    /*
        unsigned int nb_cells = ( size[0]-1 )*( size[1]-1 )*( size[2]-1 );
        cout << "valid positions: " << ( double )nb_valid/( double )nb_cells*100
            << "%" << endl;
        cout << "wrong evals: " << ( double )nb_wrong_evals/( double )nb_cells*100
            << "%" << endl;
        cout << "wrong eigenspace: "
            << ( double )nb_wrong_espace/( double )nb_cells*100
            << "%" << endl;
        cout << "wrong basis: " << ( double )nb_wrong_basis/( double )nb_cells*100
            << "%" << endl;
        cout << "percentage of candidate voxels: "
            << ( double )nbok/( double )nb_cells
            << endl;
    */
}


#endif

