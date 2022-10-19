#include "crease.hpp"
#include "face.hpp"
#include "measure_wrapper.hpp"
#include "grid.hpp"
#include "newton_search.hpp"
#include "extractor.hpp"
#include <vector>
#include <teem/hest.h>
#include <iostream>
#include <fstream>

namespace {
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

struct pvo_solution {
    pvo_solution() : p(), local_fid(0), global_fid() {}
    
    vec3 p;
    unsigned int local_fid;
    FaceId global_fid;
};

}

int main(int argc, char* argv[])
{
    using namespace spurt::crease;
    is_ridge = true;
    crease_kind = 2;
    upsample = 1;
    max_depth = 4;
    value_threshold = 0.9;
    value_threshold_select = 0.99995;
    strength_threshold = 0;
    strength_threshold_select = -1;
    
    Nrrd* nrrd = readNrrd(argv[1]);
    
    
    std::cout << "we are extracting "
              << (is_ridge ? "RIDGES" : "VALLEYS") << std::endl;
    std::cout << "crease_kind = " << crease_kind << std::endl;
    std::cout << "spatial upsampling: " << crease::upsample << std::endl;
    std::cout << "subdividion depth = " << max_depth << std::endl;
    std::cout << "threshold on value = " << value_threshold << std::endl;
    std::cout << "threshold on crease strength = " << strength_threshold << std::endl;
    
    apply_filter = true;
    nb_crossings = 0;
    
    // check if this is a tensor field <-- do we need that?
    bool is_tensor = (nrrd->dim == 4 && nrrd->axis[0].size == 7);
    
    Grid sample_grid(nrrd, crease::upsample);
    
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
    
    // create a jittered grid
#ifdef JITTER
    srand48(time(0));
    double dx = 0.1 * drand48();
    double dy = 0.1 * drand48();
    double dz = 0.1 * drand48();
#else
    double dx = 0.05; // "deterministic jittering" for reproducible results
    double dy = 0.05;
    double dz = 0.05;
#endif
    sample_grid.min[0] += dx;
    sample_grid.min[1] += dy;
    sample_grid.min[2] += dz;
    
    unsigned int nb_wrong_val = 0;
    unsigned int nb_wrong_eval = 0;
    
    // voxel data
    std::vector< vec3 > v(4); // vertex coordinates
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
    if (read_info) {
        nb_ok_vertices = read_vertex_info(ok_vertex, conf_ok);
    } else {
        nb_ok_vertices = compute_vertex_info(sample_grid, ok_vertex, conf_ok);
    }
    std::cout << "there are " << nb_ok_vertices << " valid vertices ("
              << 100*nb_ok_vertices / (M*N*P) << "%)" << std::endl;
    std::cout << std::flush;
    
    unsigned int base_eval_number = the_wrapper->get_nb_measures();
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
    boost::timer timer;
    
    for (unsigned int f = 0 ; f < n3 ; f++) {
        FaceId fid = face_id(all_faces[f].first, all_faces[f].second, M, N);
        unsigned int i, j, k, axis;
        ijk(i, j, k, all_faces[f].first, M, N);
        axis = all_faces[f].second;
        std::string label;
        if (f < n1) {
            label = "Z = cst";
        } else if (f < n2) {
            label = "Y = cst";
        } else {
            label = "X = cst";
        }
        
        int pct = 100 * f / n3;
        if (pct > prev_pct) {
            std::cout << std::endl << "progress: " << pct << "% complete - "
                      << the_wrapper->get_nb_measures() - base_eval_number << " measures so far"
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
            ids[l] = x + M * (y + N * z);
            if (!conf_ok[ids[l]]) {
                // we skip faces that contain low confidence values
                should_be_processed = false;
                break;
            } else if (ok_vertex[ids[l]]) {
                should_be_processed = true;
            }
        }
        if (!should_be_processed) {
            continue;
        }
        
        for (unsigned int l = 0 ; l < 4 ; l++) {
            ok_faces.push_back(v[l]);
        }
        
        // do the actual work
        current_face_id = fid;
        std::vector< vec3 > points;
        current_face_info.set_info(v[0], v[1], v[2], v[3]);
        prev = the_wrapper->get_nb_measures();
        bool found_something;
        
        timer.restart();
        int res = search_face(points, v[0], v[1], v[2], v[3], max_depth, found_something);
        double dt = timer.elapsed();
        
        if (res) {
            total_found_face_processing_time += dt;
            ++nb_found_face;
        } else {
            total_empty_face_processing_time += dt;
            ++nb_empty_face;
        }
        cur = the_wrapper->get_nb_measures();
        
        if (res) {
            std::cout << "extractor: PVO successful in " << dt << "sec. on face " << fid
                      << " after " << cur - prev << " measures" << std::endl;
            vec3 refpos = points.front();
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
        } else if (found_something) {
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
    std::map< FaceId, vec3 > tangents;
    
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
                } else {
                    found_invalid_solution = true;
                }
            }
        }
        
        unsigned int np = fpids.size();
        if (np == 0) {
            is_fixed[v] = true;
            continue;
        } else if (np == 1 && !found_invalid_solution) {
            std::cout << "voxel #" << v << " will be processed later" << std::endl;
            voxels_to_look_into.push_back(v);
            ++nb_bad;
            continue;
        } else if (np == 2) {
            ++nb_good;
            std::cout << "segment found in voxel #" << v << std::endl;
            vec3 p0 = all_face_points[fpids[0]];
            vec3 p1 = all_face_points[fpids[1]];
            
            // debug code begins...
            round1.push_back(p0);
            round1.push_back(p1);
            // ...debug code ends
            
            vec3 dpdt = p0 - p1;
            dpdt /= norm(dpdt);
            add_segment(fpids);
            tangents[face_ids[0]] = dpdt;
            tangents[face_ids[1]] = -1 * dpdt; // orientation points to next voxel
            
            is_fixed[v] = true;
        } else if (np > 2) {
            ++nb_good;
            // connect points with highest (resp. lowest) value
            std::vector< double > vals(np);
            for (unsigned int n = 0 ; n < np ; ++n) {
                vals[n] = the_wrapper->value(all_face_points[fpids[n]]);
            }
            std::vector< unsigned int > sorted(np);
            spurt::sort(vals, sorted);
            unsigned int id0, id1;
            if (is_ridge) {
                id0 = sorted[np-1];
                id1 = sorted[np-2];
            } else {
                id0 = sorted[0];
                id1 = sorted[1];
            }
            
            std::cout << "TRICKY segment found in voxel #" << v << std::endl;
            vec3 p0 = all_face_points[fpids[id0]];
            vec3 p1 = all_face_points[fpids[id1]];
            
            // debug code begins...
            round1.push_back(p0);
            round1.push_back(p1);
            // ...debug code ends
            
            vec3 dpdt = p0 - p1;
            dpdt /= norm(dpdt);
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
            
            cur = the_wrapper->get_nb_measures();
            ++nb_proc;
            
            unsigned int v = *it;
            unsigned int k = v / (M - 1) / (N - 1);
            unsigned int j = (v - k * (M - 1) * (N - 1)) / (M - 1);
            unsigned int i = v % (M - 1);
            
            std::vector< unsigned int > fpids;
            std::vector< unsigned int > point_face;
            std::vector< FaceId > face_ids;
            
            bool found_invalid_face = false;
            
            if (display_debug_info) {
                std::cout << "looping over faces..." << std::endl;
            }
            for (unsigned int f = 0 ; f < 6 ; f++) {
                FaceId fid = global_face_id(i, j, k, sample_grid, f);
                std::map< FaceId, unsigned int >::const_iterator it = face_point_ids.find(fid);
                if (it != face_point_ids.end()) {
                    if (it->second != invalid_face_id) {
                        if (display_debug_info) {
                            std::cout << fid << " contains pos #" << it->second << std::endl;
                        }
                        fpids.push_back(it->second);
                        point_face.push_back(f);
                        face_ids.push_back(fid);
                    } else {
                        found_invalid_face = true;
                    }
                } else {
                    if (display_debug_info) {
                        std::cout << fid << " is empty" << std::endl;
                    }
                }
            }
            
            fixing_voxel = false;
            unsigned int np = fpids.size();
            if (np == 0) {
                // that cannot happen, by construction. let the user know that s/he has switched
                // to the twilight zone
                if (display_debug_info) {
                    std::cout << "reached an empty voxel in second voxel loop!" << std::endl;
                }
                is_fixed[v] = true;
                
                the_wrapper->turn_buffer_off();
                continue;
            } else if (np == 2) {
                // this case corresponds to a secondary fix/addition of a voxel as a side effect
                // of fixing its neighbor(s): we have killed two birds with one stone!
                std::cout << "segment found in voxel #" << v << std::endl;
                vec3 p0 = all_face_points[fpids[0]];
                vec3 p1 = all_face_points[fpids[1]];
                
                // debug code begins...
                round12.push_back(p0);
                round12.push_back(p1);
                // ...debug code ends
                
                vec3 dpdt = p0 - p1;
                dpdt /= norm(dpdt);
                add_segment(fpids);
                tangents[face_ids[0]] = dpdt;
                tangents[face_ids[1]] = -1. * dpdt; // orientation points to next voxel
                is_fixed[v] = true;
            } else if (np > 2) {
                // connect points with highest (resp. lowest) value
                std::vector< double > vals(np);
                for (unsigned int n = 0 ; n < np ; ++n) {
                    vals[n] = the_wrapper->value(all_face_points[fpids[n]]);
                }
                std::vector< unsigned int > sorted(np);
                spurt::sort(vals, sorted);
                unsigned int id0, id1;
                if (is_ridge) {
                    id0 = sorted[np-1];
                    id1 = sorted[np-2];
                } else {
                    id0 = sorted[0];
                    id1 = sorted[1];
                }
                
                std::cout << "TRICKY segment found in voxel #" << v << std::endl;
                vec3 p0 = all_face_points[fpids[id0]];
                vec3 p1 = all_face_points[fpids[id1]];
                
                // debug code begins...
                round1.push_back(p0);
                round1.push_back(p1);
                // ...debug code ends
                
                vec3 dpdt = p0 - p1;
                dpdt /= norm(dpdt);
                add_segment(fpids);
                tangents[face_ids[id0]] = dpdt;
                tangents[face_ids[id1]] = -1 * dpdt; // orientation points to next voxel
                
                is_fixed[v] = true;
            } else if (np == 1) {
                if (found_invalid_face) {
                    is_fixed[v] = true;
                    continue;
                }
                
                // debug code begins...
                if (display_debug_info) {
                    std::cout << std::endl << "must now fix voxel #" << v << std::endl;
                }
                // ...debug code ends
                
                bool fixed = false;
                
                timer.restart();
                
                // build up voxel geometry
                std::vector< vec3 > vert(8);
                sample_grid.voxel(vert, i, j, k);
                if (display_debug_info)
                    std::cout << "voxel is: 0=" << vert[0] << ", 1=" << vert[1] << ", 2=" << vert[2]
                              << ", 3=" << vert[3] << ", 4=" << vert[4] << ", 5=" << vert[5] << ", 6=" << vert[6]
                              << ", 6=" << vert[6] << ", 7=" << vert[7] << std::endl;
                              
                fixing_voxel = true;
                prev = the_wrapper->get_nb_measures();
                
                std::pair< int, vec3 > first_guess; // face id / position on face
                first_guess.first = -1;   // no first guess so far
                
                // check if we have a tangent available at this position
                std::map< FaceId, vec3 >::iterator tang_it = tangents.find(face_ids[0]);
                
                // tracking solution specific to crease lines of mode
                double h = 0.25;
                if (crease_kind == 2) {
                    unsigned int fid_out;
                    vec3 out;
                    if (track_ridge_line(out, fid_out, vert, all_face_points[fpids[0]], point_face[0])) {
                        double val = the_wrapper->value(out);
                        if (display_debug_info)
                            std::cout << "tracking was successful. end position = " << out
                                      << ", face id = " << fid_out << ", value = " << val << std::endl;
                        first_guess.first = fid_out;
                        first_guess.second = out;
                        
                        // tracking approach provides us with greater accuracy
                        h = 0.1;
                    } else {
                        if (display_debug_info) {
                            std::cout << "tracking failed" << std::endl;
                        }
                    }
                } else if (tang_it != tangents.end()) {
                    unsigned int fid_out;
                    vec3 out;
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
                std::vector< vec3 > points;
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
                
                total_setup_time += timer.elapsed();
                ++nb_setup;
                
                timer.restart();
                // if we have identified a first guess, look for it in priority
                if (!found && first_guess.first >= 0) {
                    fixing_voxel = true;
                    if (display_debug_info) {
                        std::cout << "we are about to scrutinize the right face" << std::endl;
                    }
                    const unsigned int* ids = faces[first_guess.first];
                    face_type face(vert[ids[0]], vert[ids[1]], vert[ids[2]], vert[ids[3]], *spurt::crease::the_wrapper);
                    
                    if (crease_kind == 2) {
                        face_type out;
                        refine_face(out, face, first_guess.second, h);
                        if (display_debug_info)
                            std::cout << std::endl << std::endl << "we are scrutinizing face: "
                                      << out.p[0] << ", " << out.p[1] << ", " << out.p[2] << ", " << out.p[3]
                                      << std::endl;
                                      
                        unsigned int n = 0;
                        if (crease_kind < 2) {
                            n = 2;
                        }
                        for (; n < 3 && !found ; ++n) {
                        
                            if (!n) {
                                if (display_debug_info) {
                                    std::cout << "first attempt: bold approach" << std::endl;
                                }
                                spurt::crease::speedup = true;
                            } else if (n == 1) {
                                if (display_debug_info) {
                                    std::cout << "first attempt failed: cautious approach on second attempt" << std::endl;
                                }
                                spurt::crease::speedup = false;
                            } else {
                                if (display_debug_info) {
                                    std::cout << "everything we tried so far failed. using a smaller area" << std::endl;
                                }
                                spurt::crease::speedup = false;
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
                    } else { // simply look on face pointed to by tangent vector
                        if (display_debug_info)
                            std::cout << std::endl << std::endl << "we are scrutinizing face: "
                                      << face.p[0] << ", " << face.p[1] << ", " << face.p[2] << ", " << face.p[3]
                                      << std::endl;
                                      
                        spurt::crease::speedup = false;
                        current_vertices.clear();
                        current_face_info.set_info(face.p[0], face.p[1], face.p[2], face.p[3]);
                        found = search_face(points, face.p[0], face.p[1], face.p[2], face.p[3], max_depth_fix,
                                            found_something);
                                            
                        FaceId _fid = global_face_id(i, j, k, sample_grid, first_guess.first);
                        if (found) {
                            the_solution.p = points[0];
                            the_solution.local_fid = first_guess.first;
                            the_solution.global_fid = _fid;
                        } else {
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
                    spurt::crease::speedup = false;
                    spurt::crease::fixing_voxel = true;
                    for (unsigned int vf = 0 ; vf < 6 && !found && !found_something; vf++) {
                        if (vf == point_face[0] || // skip position that we already know
                                (crease_kind < 2 && // or face we have already ruled out
                                 first_guess.first >= 0 && vf == first_guess.first)) {
                            continue;
                        }
                        fixing_voxel = true;
                        FaceId fid = global_face_id(i, j, k, sample_grid, vf);
                        if (scrutinized.find(fid) != scrutinized.end()) {
                            continue;
                        }
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
                        } else if (found_something || __d == max_depth_fix) {
                            // we have reached the maximum allowed subdivision depth
                            // without being successful. we won't be any more successful
                            // in further attempts
                            scrutinized.insert(fid);
                        }
                    }
                }
                
                the_wrapper->turn_buffer_off();
                
                double dt = timer.elapsed();
                
                if (found) {
                    // add new entry to central database
                    std::cout << "found a crease point that was hiding #"
                              << added_vertices.size() + 1 << "! : " << the_solution.p << std::endl;
                    std::cout << "it took " << the_wrapper->get_nb_measures() - cur << " measures"
                              << " and a subdivision depth = " << __d << std::endl;
                    all_face_points.push_back(the_solution.p);
                    face_point_ids[the_solution.global_fid] = all_face_points.size() - 1;
                    std::cout << "we added pos #" << all_face_points.size() - 1 << " to face " << the_solution.global_fid << std::endl;
                    added_vertices.push_back(all_face_points.size() - 1);
                    std::cout << the_wrapper->get_nb_measures() << " measures so far" << std::endl;
                    fpids.push_back(all_face_points.size() - 1);
                    point_face.push_back(the_solution.local_fid);
                    face_ids.push_back(the_solution.global_fid);
                    
                    // update tangent information
                    vec3 p0 = all_face_points[fpids[0]];
                    vec3 p1 = all_face_points[fpids[1]];
                    
                    // debug code begins...
                    round2.push_back(p0);
                    round2.push_back(p1);
                    // ...debug code ends
                    
                    vec3 dpdt = p0 - p1;
                    dpdt /= norm(dpdt);
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
                } else if (found_something) {
                
                    // nothing to be done for this voxel where a solution exists, as required by
                    // the continuity of the crease lines, but does not meet our filtering criteria
                } else  {
                    if (display_debug_info) {
                        std::cout << "we failed to fix this voxel" << std::endl;
                        std::cout << "vertices are: " << std::endl;
                    }
                    std::cout << "FAILED" << std::endl;
                    for (unsigned int l = 0 ; l < 8 ; l++) {
                        problematic_voxels.push_back(vert[l]);
                        if (display_debug_info) {
                            std::cout << "#" << l << ": " << vert[l] << std::endl;
                        }
                    }
                    
                    // failed_voxels.push_back(v);
                }
                
                if (found || found_something) {
                    total_fixed_voxel_processing_time += dt;
                    ++nb_fixed_voxel;
                } else {
                    total_failed_voxel_processing_time += dt;
                    ++nb_failed_voxel;
                }
            }
        }
        // voxels_to_look_into.swap(failed_voxels);
        //  failed_voxels.clear();
    }
    std::cout << "average time spent on fixed voxels = " << total_fixed_voxel_processing_time / (float)nb_fixed_voxel
              << " sec." << std::endl
              << "average time spent on failed voxels = " << total_failed_voxel_processing_time / (float)nb_failed_voxel
              << " sec." << std::endl;
              
    std::cout << "the complete computation required " << the_wrapper->get_nb_measures() - base_eval_number
              << " measures" << std::endl;
              
              
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
              << "number of PVO computations performed: " << spurt::crease::nb_pvo << std::endl;
              
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
        
        if (a == -1) {
            connected_to[f0].first = f1;
        } else {
            connected_to[f0].second = f1;
        }
        
        if (b == -1) {
            connected_to[f1].first = f0;
        } else {
            connected_to[f1].second = f0;
        }
    }
    
    creases.clear();
    std::vector< bool > inserted(all_face_points.size(), false);
    for (unsigned int i = 0 ; i < all_edges.size() ; i++) {
        // edge end points
        unsigned int i0, i1;
        i0 = all_edges[i].first;
        i1 = all_edges[i].second;
        
        assert(i0 < all_face_points.size() && i1 < all_face_points.size());
        
        if (inserted[i0] || inserted[i1]) {
            continue;
        }
        
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
            if ((unsigned int)link0 != prev && !inserted[link0]) {
                my_list.push_back(link0);
            } else if (link1 >= 0 && (unsigned int)link1 != prev && !inserted[link1]) {
                my_list.push_back(link1);
            } else {
                break;
            }
            
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
            if ((unsigned int)link0 != prev && !inserted[link0]) {
                my_list.push_front(link0);
            } else if (link1 >= 0 && (unsigned int)link1 != prev && !inserted[link1]) {
                my_list.push_front(link1);
            } else {
                break;
            }
            
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



