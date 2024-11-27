#include <maps/symplectic4d.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <image/nrrd_wrapper.hpp>
#include <math/stat.hpp>
#include <data/locator.hpp>
#include <array>
#include <random>
#include <thread>
#include <iostream>
#include <sfcnn.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
std::atomic<size_t> tbb_progress_counter;

int nsamples, max_iter, nhits, min_hits, dim;
double thickness, k1, k2, eps, max_mean_gap;
std::string filename;

using namespace spurt;

typedef symplectic4D map_type;
typedef map_type::state_type state_type;
typedef map_type::bounds_type bounds_type;
typedef point_locator<state_type, long> NNlocator_type;
typedef NNlocator_type::point_type point_type;
typedef std::vector<state_type> orbit_type;

std::pair<double, double> compute_gaps(orbit_type& orb) {
    // ProgressDisplay timer(false);
    // timer.start();
    std::vector<point_type> points;
    for (int i=0; i<orb.size(); ++i) {
        points.push_back(point_type(orb[i], i));
    }
    NNlocator_type locator(points.begin(), points.end());
    // timer.stop();
    // std::cout << "KD tree built in " << timer.cpu_time() << " s. (cpu) / "
        // << timer.wall_time() << " s. (wall)\n";
    std::vector<point_type> nns;
    std::vector<double> gaps;
    // timer.start();
    for (auto it=orb.begin(); it!=orb.end(); ++it) {
        nns.clear();
        locator.find_n_nearest_points(nns, *it, 3);
        
        gaps.push_back(norm(*it-nns.back().position()));
    }
    // timer.stop();
    // std::cout << "Gaps of " << orb.size() << " points computed in " << timer.cpu_time()
        // << " s. (cpu) / " << timer.wall_time() << " s. (wall)\n";
    return meanvariance(gaps);
}

int main(int argc, const char* argv[]) {

    namespace cl = command_line;

    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group"),
        map_group(false, false, "Map group"),
        slab_group(false, false, "Slab group");

    cl::option_parser parser(argv[0],
        "4D Symplectic Map Visualization Using Slab Technique");

    try {
        // parser.use_default_symbol();
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("k1", k1, -2.25, "K1", map_group);
        parser.add_value("k2", k2, -3.0, "K2", map_group);
        parser.add_value("epsilon", eps, 1, "Epsilon12", map_group);
        parser.add_value("thickness", thickness, 1.0e-4, "Slab thickness", slab_group);
        parser.add_value("dimension", dim, 3, "Slice dimension", slab_group);
        parser.add_value("hits", nhits, 4000, "Number of hits", slab_group);
        parser.add_value("max", max_iter, 1000000, "Max iterations", slab_group);
        parser.add_value("min", min_hits, 200, "Min number if hits", slab_group);
        parser.add_value("gap", max_mean_gap, -1, "Max mean gap", slab_group);
        parser.add_value("samples", nsamples, 50, "Number of map samples", slab_group);
        parser.add_value("output", filename, "Output filename", required_group);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options enteredso far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }

    map_type amap(k1, k2, eps);
    slab_map slab(dim, thickness);
    typedef std::vector<state_type> hits_type;

    const bounds_type bounds = amap.bounds();
    auto min = bounds.min();
    auto max = bounds.max();
    std::vector<state_type> seeds(nsamples);

    std::random_device r;
    auto val = r();
    std::cout << "seed is " << val << '\n';
    std::mt19937 gen(val); // Standard mersenne_twister_engine seeded with r()
    std::uniform_real_distribution<> uniform(0, 1);

    for (int n=0; n<nsamples; ++n) {
        seeds[n][dim] = 0;
        for (int d=0; d<4; ++d) {
            if (d==dim) continue;
            double v = uniform(gen);
            seeds[n][d] = min[d] + v*(max[d]-min[d]);
        }
    }

    ProgressDisplay progress(true);
    progress.start(nsamples, "Sampling 4D symplectic map");
    tbb_progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<int>(0,nsamples),
                       [&](tbb::blocked_range<int> r) {

        for (int n=r.begin(); n!=r.end(); ++n) {
            std::vector<state_type> hits;
            std::vector<state_type> orbit;
            progress.update(tbb_progress_counter);
            slab.run(hits, orbit, amap, seeds[n], max_iter, nhits);
            ++tbb_progress_counter;

            if (hits.size() < min_hits) continue;

            std::pair<double, double> meanvar = compute_gaps(hits);
            const double& mean = meanvar.first;
            if (max_mean_gap > 0 && mean > max_mean_gap) continue;

            std::array<int, 2> dims{ {4, (int)hits.size()} };
            std::string name = filename + "_hits_n=" + std::to_string(n)
                + "_g=" + std::to_string(mean) + "_s=" + std::to_string(val) + ".nrrd";
            nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&hits[0]), name, dims);

            dims[1] = (int)orbit.size();
            name = filename + "_whole_orbit_n=" + std::to_string(n)
                + "_g=" + std::to_string(mean) + "_s=" + std::to_string(val) + ".nrrd";
            nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&orbit[0]), name, dims);

            std::cout << "\nhits contains " << hits.size() << " points and "
                << "orbit contains " << orbit.size() << " points\n";

        }
    });
    // for (int n=0; n<nsamples; ++n) {
    //     progress.update(n);
    //     orbit_type& orbit = all_orbits[n];
    //     slab.run(orbit, amap, seeds[n], max_iter, nhits);
    // }
    progress.stop();

    // int nhits = 0;
    // int npoints = 0;
    // for (int n=0; n<nsamples; ++n) {
    //     npoints += all_orbits[n].size();
    //     nhits += all_hits[n].size();
    // }
    //
    // double* hits = (double*)calloc(nhits*5, sizeof(double));
    // double* points = (double*)calloc(npoints*5, sizeof(double));
    // int hid, oid;
    // for (int n=0, hid=0, oid=0; n<nsamples; ++n) {
    //     const hits_type& h = all_hits[n];
    //     for (int k=0; k<h.size(); ++k, ++hid) {
    //         const state_type& s = h[k];
    //         std::cout << "s=" << s << '\n';
    //         hits[5*hid  ] = s[0];
    //         hits[5*hid+1] = s[1];
    //         hits[5*hid+2] = s[2];
    //         hits[5*hid+3] = s[3];
    //         hits[5*hid+4] = n;
    //     }
    //     const orbit_type& o = all_orbits[n];
    //     for (int k=0; k<o.size(); ++k, ++oid) {
    //         const state_type& s = o[k];
    //         points[5*oid  ] = s[0];
    //         points[5*oid+1] = s[1];
    //         points[5*oid+2] = s[2];
    //         points[5*oid+3] = s[3];
    //         points[5*oid+4] = n;
    //     }
    // }
    //
    // std::cout << "Blah" << std::endl;
    // std::array<int, 2> dims{ {5, nhits} };
    // std::string name = filename + "_hits_n=" + std::to_string(nsamples) + "_it=" + std::to_string(max_iter) + "_s=" + std::to_string(val) + ".nrrd";
    // if (nhits != 0) {
    //     nrrd_utils::writeNrrdFromContainers(hits, name, dims);
    // }
    // std::cout << "Blih" << std::endl;
    // dims[1] = npoints;
    // name = filename + "_orbits_n=" + std::to_string(nsamples) + "_it=" + std::to_string(max_iter) + "_s=" + std::to_string(val) + ".nrrd";
    // nrrd_utils::writeNrrdFromContainers(points, name, dims);
    return 0;
}
