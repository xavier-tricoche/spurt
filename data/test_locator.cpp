#include "locator.hpp"
#include <iostream>
#include <string>
#include <misc/time_helper.hpp>

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
bool test_solution(const std::vector<P>& allpoints, const P& query, const P& answer)
{
    double d = norm(answer.coordinate()-query.coordinate());
    double secondbest = std::numeric_limits<double>::max();
    int secondbesti = -1;
    for (int i=0 ; i<allpoints.size() ; ++i) {
        if (i == answer.data()) {
            continue;
        }
        double dd = norm(allpoints[i].coordinate()-query.coordinate());
        if (dd < d) {
            return false;
        } else if (dd < secondbest) {
            secondbest = dd;
            secondbesti = i;
        }
    }
    if (drand48() < 0.01) {
        std::cout << "closest point to " << query.coordinate() << " is "
                  << answer.coordinate() << " at distance " << d << " (" << 100*d/sqrt(8) << "%)" << std::endl;
        std::cout << "second closest point was " << allpoints[secondbesti].coordinate()
                  << " at distance " << secondbest << std::endl;
    }
    return true;
}

template<int K>
void test_locator(int n, int ns)
{
    typedef spurt::point_locator<double, int, K>   locator_type;
    typedef typename locator_type::point_type       point_type;
    typedef typename locator_type::coord_type       coord_type;
    
    std::vector<point_type> all_points(n);
    locator_type locator;
    srand48(time(0));
    for (int i=0 ; i<n ; ++i) {
        coord_type c;
        for (int k=0 ; k<K ; ++k) {
            c[k] = -1. + drand48() * 2.;
        }
        all_points[i] = point_type(c, i);
        locator.insert(all_points[i]);
    }
    
    spurt::timer _timer;
    double total_time = 0;
    int nbwrong = 0;
    for (int i=0 ; i<ns ; ++i) {
        coord_type c;
        for (int k=0 ; k<K ; ++k) {
            c[k] = -1. + drand48() * 2.;
        }
        _timer.start();
        point_type nearest = locator.find_nearest_point(c);
        total_time += _timer.elapsed();
        if (!test_solution(all_points, point_type(c, 0), nearest)) {
            ++nbwrong;
        }
    }
    
    std::cout << "Results:\n"
              << "Number of incorrect answers: " << nbwrong << " (" << 100*nbwrong/ns << "%)\n"
              << "Average time per query: " << total_time << " s.\n";
    exit(0);
}

int main(int argc, char* argv[])
{
    int n = 1000;
    int d = 3;
    int ns = 10000;
    for (int i=1 ; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-n" || arg == "--number") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing number");
            }
            n = atoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            printUsageAndExit(argv[0]);
        } else if (arg == "-d" || arg == "--dimension") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing dimension");
            }
            d = atoi(argv[++i]);
        } else if (arg == "-ns" || arg == "--nsamples") {
            if (i == argc-1) {
                printUsageAndExit(argv[0], "missing number of samples");
            }
            ns = atoi(argv[++i]);
        }
    }
    
    switch (d) {
        case 1:
            test_locator<1>(n, ns);
        case 2:
            test_locator<2>(n, ns);
        case 3:
            test_locator<3>(n, ns);
        case 4:
            test_locator<4>(n, ns);
        case 5:
            test_locator<5>(n, ns);
        case 6:
            test_locator<6>(n, ns);
        default:
            std::cerr << "this test program only covers dimensions 1 through 6\n";
    }
    return 0;
}
