#include "image/creaseline2d.hpp"

std::vector< double > xavier::crease::crease_strength;
std::vector< unsigned int > xavier::crease::isolated;
std::vector< bool > xavier::crease::skipped;
std::vector< bool > xavier::crease::singular;
float* xavier::crease::reconstructed_image = 0;

struct PosOrder {
    bool operator()(const nvis::vec2& p0, const nvis::vec2& p1) const {
        return (p0[0] > p1[0] + 0.00001 ||
                (fabs(p0[0] - p1[0]) < 0.00001 && p0[1] > p1[1] + 0.00001));
    }
};

void uniquify(std::map< unsigned int, std::pair< int, int > >& cell2points,
              std::vector< nvis::vec2 >& pos)
{
    unsigned int N = pos.size();
    std::map< nvis::vec2, unsigned int, PosOrder > ref;
    std::map< nvis::vec2, unsigned int, PosOrder >::iterator it;
    std::vector< unsigned int > old2new(pos.size());
    std::vector< nvis::vec2 > unique;
    std::vector< double > str;
    unsigned int id = 0;
    for (unsigned int i = 0 ; i < pos.size() ; ++i) {
        it = ref.find(pos[i]);
        if (it != ref.end()) {
            old2new[i] = it->second;
        }
        else {
            old2new[i] = id;
            ref[pos[i]] = id;
            unique.push_back(pos[i]);
            str.push_back(xavier::crease::crease_strength[i]);
            ++id;
        }
    }

    for (unsigned int i = 0 ; i < cell2points.size() ; ++i) {
        int i0 = cell2points[i].first;
        int i1 = cell2points[i].second;
        if (i0 >= 0) {
            cell2points[i].first = old2new[i0];
        }
        if (i1 >= 0) {
            cell2points[i].second = old2new[i1];
        }
    }

    std::swap(pos, unique);
    std::swap(xavier::crease::crease_strength, str);

    std::cout << "after uniquification: " << pos.size() << " positions instead of "
              << N << std::endl;
}

void connect_segments(std::vector< std::list< unsigned int > >& creases,
                      const std::map< unsigned int, std::pair< int, int > >& cellpoints,
                      unsigned int npoints)
{
    // compute point <- point -> point connections
    std::vector< std::pair< int, int > > connected_to(npoints, std::pair< int, int >(-1, -1));

    std::map< unsigned int, std::pair< int, int > >::const_iterator mit;
    for (mit = cellpoints.begin() ; mit != cellpoints.end() ; ++mit) {

        int id0 = mit->second.first;
        int id1 = mit->second.second;

        if (id0 >= 0 && id1 >= 0) {
            if (connected_to[id0].first < 0) {
                connected_to[id0].first = id1;
            }
            else {
                connected_to[id0].second = id1;
            }

            if (connected_to[id1].first < 0) {
                connected_to[id1].first = id0;
            }
            else {
                connected_to[id1].second = id0;
            }
        }
    }

    creases.clear();
    std::vector< bool > inserted(npoints, false);
    for (mit = cellpoints.begin() ; mit != cellpoints.end() ; ++mit) {

        // edge end points
        int i0, i1;
        i0 = mit->second.first;
        i1 = mit->second.second;
        if (i0 < 0 || i1 < 0) continue;

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
            assert(link1 < 0 || link1 < npoints);
            if ((unsigned int)link0 != prev && !inserted[link0])
                my_list.push_front(link0);
            else if (link1 >= 0 && (unsigned int)link1 != prev && !inserted[link1])
                my_list.push_front(link1);
            else break;

            prev = cur;
            cur = my_list.front();
        }

        creases.push_back(my_list);
    }

    std::cout << "verifying results:" << std::endl;
    for (unsigned int i = 0 ; i < creases.size() ; i++) {
        std::cout << "component #" << i << ": (" << creases[i].size() << ")" << std::flush;
        std::list< unsigned int >::iterator lit;
        for (lit = creases[i].begin() ; lit != creases[i].end() ; ++lit) {
            std::cout << *lit << " " << std::flush;
        }
        std::cout << std::endl;
    }
}

void xavier::crease::
extract(const Nrrd* nrrd,
        const raster_grid<2>& grid,
        double threshold,
        bool ridge,
        std::vector< nvis::vec2 >& intersections,
        std::vector< std::list< unsigned int > >& creases)
{
    using namespace nvis;
    using namespace datastructure;
    using namespace gage_interface;

    skipped.resize(grid.resolution()[0]*grid.resolution()[1]);
    std::fill(skipped.begin(), skipped.end(), false);

    singular.resize(grid.resolution()[0]*grid.resolution()[1]);
    std::fill(singular.begin(), singular.end(), false);

    if (reconstructed_image)
        delete[] reconstructed_image;
    reconstructed_image = (float *)calloc(grid.resolution()[0] * grid.resolution()[1], sizeof(float));

    crease_strength.clear();
    isolated.clear();

    double val0, val1;
    std::vector< vec2 > x(5);
    std::vector< vec2 > p(5);
    vec2 inter;

    crease::nsub = 0;
    scalar_wrapper gH(nrrd);

    std::map< unsigned int, std::pair< int, int > > cellpoints; // crease points / cell
    std::map< Edge, unsigned int > edge2point; // crease point / edge
    std::vector< unsigned int > allpointsincell;
    std::set< Edge > processed;
    Edge k[4];

    intersections.clear();
    creases.clear();

    std::cout << "input grid has " << grid.resolution()[0] << " x " << grid.resolution()[1] << " cells\n";
    std::cout << "we are extracting RIDGES\n";

    // loop over all cells
    int pct = -1;
    for (unsigned int j = 0 ; j < grid.resolution()[1] - 1 ; j++)
        for (unsigned int i = 0 ; i < grid.resolution()[0] - 1 ; i++) {

            unsigned int n = i + j * grid.resolution()[0];
            int new_pct = floor(100 * n / grid.resolution()[0] / grid.resolution()[1]);
            if (new_pct > pct) {
                std::cout << (pct < 10 ? "0" : "") << new_pct << "\% completed\r" << std::flush;
                pct = new_pct;
            }

            // cyclic list of vertices
            p[0] = grid(i + 1, j);
            p[1] = grid(i + 1, j + 1);
            p[2] = grid(i, j + 1);
            p[3] = grid(i, j);
            p[4] = p[0];

            nvis::vec2 __p = 0.25 * (p[0] + p[1] + p[2] + p[3]);
            // __p[0] *= (double)(grid.resolution()[0] - 1);
            // __p[1] *= (double)(grid.resolution()[1] - 1);
            double __v;
            gH.value(__p, __v);
            reconstructed_image[i+j*grid.resolution()[0]] = __v;
            // std::cout << "image(" << i << ", " << j << ") = " << __v << std::endl;


            // corresponding cyclic list of edges as a [ cell <-- edge --> cell ]
            // relationship
            k[0] = Edge(cellid(i, j, grid), right(i, j, grid));
            k[1] = Edge(cellid(i, j, grid), top(i, j, grid));
            k[2] = Edge(cellid(i, j, grid), left(i, j, grid));
            k[3] = Edge(cellid(i, j, grid), down(i, j, grid));

            // loop over 4 edges
            allpointsincell.clear();
            for (unsigned int l = 0 ; l < 4 ; l++) {
                std::map< Edge, unsigned int >::iterator _it = edge2point.find(k[l]);
                if (_it != edge2point.end()) {
                    allpointsincell.push_back(_it->second);
                    // std::cout << "reusing crease point #" << allpointsincell.back()
                    // << std::endl;
                    continue;
                }
                else if (processed.find(k[l]) != processed.end()) {
                    // std::cout << "(" << k[l]._i << ", " << k[l]._j << ") processed already - skipping\n";
                    continue;
                }
                processed.insert(k[l]);

                // Nrrd raster is scaled to [0,resolution()[0]-1] x [0,resolution()[1]-1]
                // that assumes that grid is invariably a unit square
                // x[l][0] = (grid.resolution()[0] - 1) * p[l][0];
                // x[l][1] = (grid.resolution()[1] - 1) * p[l][1];
                x[l] = p[l];
                if (!gH.value(x[l], val0)) {
                    // std::cout << "wrong value at " << x[l] << "\n";
                    skipped[i+j*grid.resolution()[0]] = true;
                    continue;
                }
                // x[l+1][0] = (grid.resolution()[0] - 1) * p[l+1][0];
                // x[l+1][1] = (grid.resolution()[1] - 1) * p[l+1][1];
                x[l+1] = p[l+1];
                if (!gH.value(x[l+1], val1)) {
                    // std::cout << "wrong value #1\n";
                    skipped[i+j*grid.resolution()[0]] = true;
                    continue;
                }
                if ((!ridge && (val0 >= threshold || val1 >= threshold)) ||
                    (ridge && (val0 <= threshold || val1 <= threshold))) {
                    skipped[i+j*grid.resolution()[0]] = true;

                    // std::cout << "skipping edge (" << val0 << ", " << val1 << ")\n";
                    continue;
                }
                double strength;
                if (find_intersection(x[l], x[l+1], inter,
                                      gH, ridge, strength)) {

                    // std::cout << "found an intersection!\n";
                    intersections.push_back(inter);
                    crease_strength.push_back(strength);
                    unsigned int interid = intersections.size() - 1;
                    edge2point[k[l]] = interid;
                    allpointsincell.push_back(interid);
                }
                else {
                    // std::cout << "no intesection on that edge\n";
                }
            }

            if (allpointsincell.size() == 1) {
                // std::cout << "found only one crease point" << std::endl;
                isolated.push_back(allpointsincell[0]);
                singular[i+j*grid.resolution()[0]] = true;
            }
            else if (allpointsincell.size() > 1) {
                unsigned int id0, id1;
                if (allpointsincell.size() == 2) {
                    id0 = allpointsincell[0];
                    id1 = allpointsincell[1];
                }
                else {
                    // std::cout << "picking best crease points out of "
                    // << allpointsincell.size() << ": ";

                    // connect points corresponding to largest / smallest values
                    // of considered scalar field
                    std::vector< double > vals(allpointsincell.size());
                    for (unsigned int k = 0 ; k < allpointsincell.size() ; k++) {
                        unsigned int pid = allpointsincell[k];
                        nvis::vec2 p = intersections[pid];
                        double val;
                        gH.value(p, val);
                        vals[k] = val;
                    }
                    double max1 = vals[0];
                    unsigned int max1i = 0;
                    double max2 = vals[1];
                    unsigned int max2i = 1;
                    for (unsigned int k = 2 ; k < allpointsincell.size() ; k++) {
                        if ((ridge && vals[k] > max1) ||
                            (!ridge && vals[k] < max1)) {
                            max1 = vals[k];
                            max1i = k;
                        }
                        else if ((ridge && vals[k] < max2) &&
                                 (!ridge && vals[k] < max2)) {
                            max2 = vals[k];
                            max2i = k;
                        }
                    }

                    // std::cout << max1i << " and " << max2i << std::endl;
                    id0 = allpointsincell[max1i];
                    id1 = allpointsincell[max2i];
                }

                cellpoints[cellid(i,j,grid)] = std::pair< int, int >(id0, id1);
            }
        }

    /*
        // we have:
        // cellpoints: cell --> pair of crease points
        // we want inverse mapping:
        // vertexcells: crease point --> pair of cells
        std::map< unsigned int, std::pair< int, int > > vertexcells;
        for (std::map< unsigned int, std::pair< int, int > >::iterator
             it = cellpoints.begin() ; it != cellpoints.end() ; it++) {
            int c = it->first; // start cell

            // corresponding crease points
            int v1, v2;
            v1 = it->second.first;
            v2 = it->second.second;
            if (v1 >= 0)
                if (vertexcells.find(v1) != vertexcells.end())
                    vertexcells[v1].second = c;
                else
                    vertexcells[v1] = std::pair< int, int >(c, -1);
            if (v2 >= 0)
                if (vertexcells.find(v2) != vertexcells.end())
                    vertexcells[v2].second = c;
                else
                    vertexcells[v2] = std::pair< int, int >(c, -1);
        }

        std::vector< bool > checked(intersections.size(), false);
        for (std::map< unsigned int, std::pair< int, int > >::iterator
             it = cellpoints.begin() ; it != cellpoints.end() ; it++) {
            int ccur, v1, v2;
            ccur = it->first; // current cell
            v1 = it->second.first; // current "left" vertex
            v2 = it->second.second; // current "right" vertex

            std::list< unsigned int > newlist;
            std::map< unsigned int, std::pair< int, int > >::iterator vit, cit;

            // current cell vertices have been previously included in a valley
            if (v1 < 0 || checked[v1] || v2 < 0) continue;

            // otherwise loop in both directions
            newlist.push_front(v1);
            while (!checked[v1]) {
                checked[v1] = true;

                vit = vertexcells.find(v1);
                if (vit == vertexcells.end()) break;

                int cnext = next(vit->second, ccur);
                if (cnext < 0) break;

                cit = cellpoints.find(cnext);
                if (cit == cellpoints.end()) break;
                ccur = cnext;

                int vnext = next(cit->second, v1);
                if (vnext < 0) break;

                newlist.push_front(vnext);
                v1 = vnext;
            }

            ccur = it->first;
            newlist.push_back(v2);
            while (!checked[v2]) {
                checked[v2] = true;

                vit = vertexcells.find(v2);
                if (vit == vertexcells.end()) break;

                int cnext = next(vit->second, ccur);
                if (cnext < 0) break;

                cit = cellpoints.find(cnext);
                if (cit == cellpoints.end()) break;
                ccur = cnext;

                int vnext = next(cit->second, v2);
                if (vnext < 0) break;

                newlist.push_back(vnext);
                v2 = vnext;
            }

            creases.push_back(newlist);
        }
    */
    uniquify(cellpoints, intersections);

    creases.clear();
    connect_segments(creases, cellpoints, intersections.size());
}

































