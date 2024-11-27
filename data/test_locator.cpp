#include "locator.hpp"
#include <iostream>
#include <string>
#include <misc/progress.hpp>
#include <misc/cxxopts.hpp>
#include <data/heap.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
std::atomic<size_t> tbb_progress_counter;
std::atomic<size_t> tbb_nb_wrong;
std::atomic<size_t> tbb_nb_tested;
std::atomic<size_t> tbb_total_time;

using namespace spurt;

void printUsageAndExit( const std::string& argv0, const std::string& offending="",
                        bool doExit = true )
{
    if (offending != "") {
        std::cerr << "ERROR: " << offending << std::endl;
    }
    std::cerr
            << "Usage  : " << argv0 << " [options]\n"
            << "Options:\n"
            << "    -h  | --help                      Print this information\n"
            << "    -n  | --number <float>            Number of points used in test\n"
            << "    -d  | --dimension <int>           Spatial dimensions used in test\n"
            << "    -ns | --nsamples <int>            Number of samples used in test\n"
            << std::endl;
            
    if (doExit) {
        exit(1);
    }
}

template<typename P>
bool test_solution(const std::vector<P>& allpoints, 
                   const P& query, const P& answer)
{
    typedef typename P::pos_type pos_type;
    typedef spurt::data_traits<pos_type> traits_type;
    double d = traits_type::norm(answer.position()-query.position());
    double secondbest = std::numeric_limits<double>::max();
    int secondbesti = -1;
    for (int i=0 ; i<allpoints.size() ; ++i) {
        if (i == answer.data()) {
            continue;
        }
        double dd = traits_type::norm(allpoints[i].position()-query.position());
        if (dd < d) {
            return false;
        } else if (dd < secondbest) {
            secondbest = dd;
            secondbesti = i;
        }
    }
    return true;
}

template<typename P>
struct further {
    typedef P point_type;
    typedef typename P::pos_type pos_type;
    pos_type m_from;
    
    further(const pos_type& from=pos_type(0)) : m_from(from) {}
    further(const point_type& from): m_from(from.position()) {}
        
    bool operator()(const point_type& p0, const point_type& p1) {
        return norm(p0.position()-m_from) > norm(p1.position()-m_from);
    }
};


template<typename P>
bool test_knn_solution(const std::vector<P>& allpoints, 
                       const P& query, size_t knn, const std::vector<P>& answer)
{
    typedef typename P::pos_type pos_type;
    // using a max distance heap with fixed capacity
    // once capacity is reached, each new push will remove current top 
    // if larger than new value
    typedef spurt::bounded_heap<P, further<P>> heap_type;
    heap_type _heap(knn, further(query));
    
    for (size_t i=0 ; i<allpoints.size(); ++i) 
    {
        _heap.push(allpoints[i], false);
    }
    
    assert(_heap.size() == knn);
    
    double maxdist = norm(_heap.top().position()-query.position());
    for (auto it=answer.begin(); it!=answer.end(); ++it)
    {
        const P& p = it->position();
        if (norm(p.position()-query.position()) > maxdist) {
            std::cout << "rejecting " << knn << "-nearest neighbor " << p.position() << " with distance " << norm(p.position()-query.position()) << " (>" << maxdist << ")\n";
            std::cout << "The contents of the heap was:\n";
            for (auto jit=_heap.begin(); jit!=_heap.end(); ++jit) {
                std::cout << "point #" << jit->data() << " " << jit->position() << " is at distance " << norm(jit->position()-query.position()) << '\n';
            }
            return false;
        }
    }
    return true;
}

template<int K>
void test_nearest_locator(int n, int ns, float f)
{
    typedef spurt::small_vector<double, K>         pos_type;
    typedef spurt::point_locator<pos_type, int>    locator_type;
    typedef typename locator_type::point_type      point_type;
    
    std::vector<point_type> all_points(n);
    srand48(time(0));
    for (int i=0 ; i<n ; ++i) {
        pos_type c;
        for (int k=0 ; k<K ; ++k) {
            c[k] = -1. + drand48() * 2.;
        }
        all_points[i] = point_type(c, i);
    }
    locator_type locator(all_points.begin(), all_points.end());
    
    tbb_total_time = 0;
    tbb_progress_counter = 0;
    tbb_nb_wrong = 0;
    tbb_nb_tested = 0;
    
    spurt::timer ttimer;
    tbb::parallel_for(tbb::blocked_range<int>(0,ns),
                       [&](tbb::blocked_range<int> r) {
        for (int p=r.begin(); p!=r.end(); ++p) {
            pos_type c;
            for (int k=0 ; k<K ; ++k) {
                c[k] = -1. + drand48() * 2.;
            }
    
            spurt::timer _timer;
            _timer.start();
            point_type nearest = locator.find_nearest_point(c);
            tbb_total_time += _timer.elapsed();
            if (drand48() < f) {
                ++tbb_nb_tested;
                if (!test_solution(all_points, point_type(c, 0), nearest)) {
                    ++tbb_nb_wrong;
                }
            }
        }
    });
    ttimer.stop();
    
    std::cout << "Results of NN test:\n"
              << "Total wall time: " << ttimer.wall_time() << " s.\n"
              << "Total CPU time: " << ttimer.cpu_time() << " s.\n"
              << "Number of incorrect answers out of " << tbb_nb_tested << " tests: " << tbb_nb_wrong << " (" << 100*float(tbb_nb_wrong)/float(tbb_nb_tested) << "%)\n"
              << "Average time per query after " << ns << " queries: " << float(tbb_total_time)/float(ns) << " s.\n";
}

template<int K>
void test_knn_locator(int n, int ns, int knn, float f)
{
    typedef spurt::small_vector<double, K>         pos_type;
    typedef spurt::point_locator<pos_type, int>    locator_type;
    typedef typename locator_type::point_type      point_type;
    
    std::vector<point_type> all_points(n);
    srand48(time(0));
    for (int i=0 ; i<n ; ++i) {
        pos_type c;
        for (int k=0 ; k<K ; ++k) {
            c[k] = -1. + drand48() * 2.;
        }
        all_points[i] = point_type(c, i);
    }
    locator_type locator(all_points.begin(), all_points.end());
    
    tbb_total_time = 0;
    tbb_progress_counter = 0;
    tbb_nb_wrong = 0;
    tbb_nb_tested = 0;
    
    spurt::timer ttimer;
    tbb::parallel_for(tbb::blocked_range<int>(0,ns),
                       [&](tbb::blocked_range<int> r) {
        for (int p=r.begin(); p!=r.end(); ++p) {
            pos_type c;
            for (int k=0 ; k<K ; ++k) {
                c[k] = -1. + drand48() * 2.;
            }
    
            spurt::timer _timer(true);
            std::vector<point_type> nns;
            locator.find_n_nearest_points(nns, c, knn);
            _timer.stop();
            tbb_total_time += _timer.cpu_time();
            if (drand48() < f) {
                ++tbb_nb_tested;
                if (!test_knn_solution(all_points, point_type(c, 0), knn, nns)) {
                    ++tbb_nb_wrong;
                }
            }
        }
    });
    ttimer.stop();
    
    std::cout << "Results of k-NN test:\n"
              << "total wall time: " << ttimer.wall_time() << " s.\n"
              << "total CPU time: " << ttimer.cpu_time() << " s.\n"
              << "Number of incorrect answers out of " << tbb_nb_tested << ": " << tbb_nb_wrong << " (" << 100*float(tbb_nb_wrong)/float(tbb_nb_tested) << "%)\n"
              << "Average time per query after " << ns << " queries: " << float(tbb_total_time)/float(ns) << " s.\n";
    exit(0);
}

int main(int argc, const char* argv[])
{   
    cxxopts::Options options("test_locator", "Test point locator for nearest neighbor and k-nearest neighbor searches");
    options.add_options()
        ("n,number", "Number of points in locator", cxxopts::value<int>()->default_value("1000"))
        ("s,samples", "Number of sample points in test", cxxopts::value<int>()->default_value("10000"))
        ("d,dimension", "Dimension of search space", cxxopts::value<int>()->default_value("3"))
        ("k", "Number of nearest neighbors to search for", cxxopts::value<int>()->default_value("8"))
        ("f,frequency", "Frequency of correctness checks", cxxopts::value<float>()->default_value("0.01"))
        ("h,help", "Print usage information");
    
    auto result = options.parse(argc, argv);
    
    if (result.count("help")) {
        std::cout << options.help() << '\n';
        exit(0);
    }
    
    int n = result["number"].as<int>();
    int ns = result["samples"].as<int>();
    int d = result["dimension"].as<int>();
    int k = result["k"].as<int>();
    float f = result["frequency"].as<float>();
    
    switch (d) {
        case 1:
            test_nearest_locator<1>(n, ns, f);
            test_knn_locator<1>(n, ns, k, f);
        case 2:
            test_nearest_locator<2>(n, ns, f);
            test_knn_locator<2>(n, ns, k, f);
        case 3:
            test_nearest_locator<3>(n, ns, f);
            test_knn_locator<3>(n, ns, k, f);
        case 4:
            test_nearest_locator<4>(n, ns, f);
            test_knn_locator<4>(n, ns, k, f);
        case 5:
            test_nearest_locator<5>(n, ns, f);
            test_knn_locator<5>(n, ns, k, f);
        case 6:
            test_nearest_locator<6>(n, ns, f);
            test_knn_locator<6>(n, ns, k, f);
        case 7:
            test_nearest_locator<7>(n, ns, f);
            test_knn_locator<7>(n, ns, k, f);
        case 8:
            test_nearest_locator<8>(n, ns, f);
            test_knn_locator<8>(n, ns, k, f);
        case 9:
            test_nearest_locator<9>(n, ns, f);
            test_knn_locator<9>(n, ns, k, f);
        case 10:
            test_nearest_locator<10>(n, ns, f);
            test_knn_locator<10>(n, ns, k, f);
        case 11:
            test_nearest_locator<11>(n, ns, f);
            test_knn_locator<11>(n, ns, k, f);
        case 12:
            test_nearest_locator<12>(n, ns, f);
            test_knn_locator<12>(n, ns, k, f);
        default:
            std::cerr << "this test program only covers dimensions 1 through 12\n";
    }
    return 0;
}