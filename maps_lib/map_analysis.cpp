#include "map_analysis.hpp"
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <data/grid.hpp>
#include <data/raster_data.hpp>

#if _OPENMP
#include <omp.h>
#endif

using namespace spurt;

typedef grid<double, 2>                         grid_type;
typedef std::pair<nvis::vec2, int>              data_type;
typedef raster_data<data_type, double, 2>       dataset_type;


double spurt::average_distance(const std::vector<nvis::vec2>& steps, int p,
                                const map_metric& metric)
{
    if (steps.size() < p) {
        return std::numeric_limits<double>::max();
    }
    double d = 0;
    int count = 0;
    for (int i=0 ; i<steps.size() ; ++i) {
        if (i+p<steps.size()) {
            d += metric.distance(steps[i], steps[i+p]);
            ++count;
        }
        if (i>=p) {
            d += metric.distance(steps[i-p], steps[i]);
            ++count;
        }
    }
    return d/(double)count;
}

int spurt::best_period(const std::vector<nvis::vec2>& steps, int maxp,
                        const map_metric& metric)
{
    std::map<double, int> period_to_dist;
    for (int p=1 ; p<=maxp ; ++p) {
        period_to_dist[average_distance(steps, p, metric)] = p;
    }
    return period_to_dist.begin()->second;
}
