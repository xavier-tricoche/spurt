#ifndef __XAVIER_PERIOD_HPP__
#define __XAVIER_PERIOD_HPP__

#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <math/fixed_vector.hpp>
#include <maps_lib/definitions.hpp>
#include <misc/sort.hpp>

#ifdef FFTW3
#include <fftw3.h>
#endif

namespace {
nvis::bbox2 bounds(const std::vector<nvis::vec2>& p)
{
    if (p.empty()) {
        return nvis::bbox2();
    }
    nvis::bbox2 bbox(p[0], p[0]);
    for (int i = 1 ; i < p.size() ; ++i) {
        bbox.add(p[i]);
    }
    return bbox;
}

const double max_angle=90;

inline bool is_full(const std::pair<int, int>& p)
{
    return (p.first >= 0 && p.second >= 0);
}

void closest(std::vector<int>& r, const std::vector<nvis::vec2>& pts, int id, int n,
             const spurt::default_metric_type& metric)
{
    r.resize(n);
    std::map<double, int> dist;
    for (int i=0 ; i<pts.size() ; ++i) {
        if (i == id) {
            continue;
        }
        dist.insert(std::pair<double, int>(metric.distance(pts[id], pts[i]), i));
    }
    int count = 0;
    for (std::map<double, int>::const_iterator i=dist.begin() ;
            i!=dist.end() && count<n ; ++i) {
        r[count++] = i->second;
    }
}

nvis::vec2 tangent(const std::vector<int>& ids, const std::vector<nvis::vec2>& pts,
                   int start, int n, const spurt::default_metric_type& metric)
{
    // compute covariance matrix
    nvis::mat2 covar(0);

    int N = std::min<int>(ids.size(), n);
    for (int i=0 ; i<N ; ++i) {
        nvis::vec2 d = metric.displacement(pts[start], pts[ids[i]]);
        d /= nvis::norm(d);
        covar += nvis::outer(d, d);
    }
    covar /= (double)N;

    // compute deviatoric
    double tr = nvis::trace(covar);
    double alpha = covar(0,0) - 0.5*tr;
    double beta = covar(0,1);
    double theta = atan((sqrt(alpha*alpha + beta*beta)-alpha)/beta);
    return nvis::vec2(cos(theta), sin(theta));
}

} // anonymous

namespace map_analysis {
#ifdef FFTW3
double fft(const std::vector<double>& q)
{
    size_t N = q.size();
    fftw_complex* in, *out;
    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for (int i = 0 ; i < N ; ++i) {
        in[i][0] = q[i];
        in[i][1] = 0.;
    }
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_free(in);

    double res = out[0][0];
    // std::cout << "fft returns " << res << std::endl;
    fftw_free(out);

    return res / (double)N;
}
#endif

inline void output_broken_chain(const std::vector<nvis::vec2>& steps,
                                const metric_type& metric)
{
    size_t size = steps.size();
    const nvis::vec2& x0 = steps.front();
    std::cerr << "\n\ndisplaying pathologic orbit seeded at " << x0 << std::endl;
    std::cerr << "orbit contains " << size << " points\n";
    std::vector< double > q(size - 1);
    std::vector< double > d(size - 1);
    for (unsigned int i = 1 ; i < size ; ++i) {
        nvis::vec2 dist = steps[i] - x0;
        q[i-1] = (double)i / dist[0] * metric.width();
        d[i-1] = metric.distance(x0, steps[i]);
        std::cerr << "step #" << i << "/" << size << ": x=" << steps[i] << ", q_approx = " << q[i] << ", d = " << d[i] << std::endl;
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));
    std::cerr << "min d = " << min << std::endl;
    double total_q = 0;
    double total_w = 0;
    for (int i = 0 ; i < q.size() ; ++i) {
        double w = min / fabs(d[i]);
        total_w += w;
        total_q += w * q[i];
        std::cerr << "weight #" << i << " = " << w << std::endl;
    }
    std::cerr << "returned value = " << total_q << "/" <<  total_w << " = "
              << total_q / total_w << std::endl;
}

inline double
min_dist_period_x_periodic(unsigned int p, const std::vector<nvis::vec2>& steps,
                           const metric_type& metric)
{
    if (steps.size() <= p) {
        return -1;
    }
    std::vector<double> pdist;
    for (int i = 0 ; i < steps.size() ; ++i) {
        if (i < p && i + p < steps.size()) {
            pdist[i] = metric.distance(steps[i], steps[i+p]);
        } else if (i >= p && i + p >= steps.size()) {
            pdist[i] = metric.distance(steps[i-p], steps[i]);
        } else if (i >= p && i + p < steps.size()) {
            pdist[i] = 0.5 * (metric.distance(steps[i], steps[i+p]) + metric.distance(steps[i-p], steps[i]));
        } else {
            pdist[i] = std::numeric_limits<double>::max();
        }
    }

    int mini = std::distance(pdist.begin(), std::max_element(pdist.begin(), pdist.end()));
    double disti = (mini + p < steps.size() ? (steps[mini+p][0] - steps[mini][0]) : (steps[mini][0] - steps[mini-p][0]));
    return (double)p / fabs(disti);
}

inline double
dist_based_x_period(const std::vector<nvis::vec2>& steps, const metric_type& metric)
{
    // NB: this method has quadratic complexity...

    typedef std::pair<unsigned int, unsigned int>   edge_id_type;
    if (steps.size() <= 1) {
        return -1;
    }

    double mind = std::numeric_limits<double>::max();
    edge_id_type min_edge_id;
    for (unsigned int i = 0 ; i < steps.size() - 1 ; ++i) {
        for (unsigned int j = i + 1 ; j < steps.size() ; ++j) {
            double d = metric.distance(steps[i], steps[j]);
            if (d < mind) {
                mind = d;
                min_edge_id = edge_id_type(i, j);
            }
        }
    }

    unsigned int i = min_edge_id.first;
    unsigned int j = min_edge_id.second;
    double q = (double)(j - i) * metric.width() / (steps[j][0] - steps[i][0]);
    // if (q < 1) {
    //  std::cerr << "i=" << i << ", j=" << j << ", dx="
    //  << steps[j][0] - steps[i][0] << ", q=" << q << '\n';
    // }
    return q;
}

std::pair<double, double>
inline period_x_periodic(const std::vector<nvis::vec2>& steps, const metric_type& metric)
{
    if (steps.size() <= 1) {
        return std::make_pair(-1, 0);
    }

    const nvis::vec2& x0 = steps.front();

    std::vector< double > q(steps.size() - 1);
    std::vector< double > d(steps.size() - 1);
    for (int i = 1 ; i < steps.size() ; ++i) {
        nvis::vec2 dist = steps[i] - x0;
        q[i-1] = (double)i / dist[0] * metric.width();
        d[i-1] = metric.distance(x0, steps[i]);
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));

    double total_q = 0;
    double total_w = 0;
    for (int i = 0 ; i < q.size() ; ++i) {
        double w = min / fabs(d[i]);
        total_w += w;
        total_q += w * q[i];
    }

    // std::cerr << "estimated period is " << total_q / total_w << '\n';
    if (std::isnan(total_q / total_w)) {
        output_broken_chain(steps, metric);
    }

    if (total_q / total_w < 0) {
        // std::cerr << "invalid period: " << total_q << " / " << total_w << std::endl;
        // std::cerr << "setting to 0\n";
        // if (x0[1]>40)
        nvis::bbox2 box = bounds(steps);
        if (true || box.max()[1] > 32) {
            output_broken_chain(steps, metric);
        }
        total_q = 0.;
    }

    return std::make_pair(total_q / total_w, 1.);
}

std::pair<double, double>
inline period_x_periodic(const std::vector<nvis::vec2>& steps,
                         const nvis::vec2& x0, const metric_type& metric)
{
    if (!steps.size()) {
        return std::make_pair(-1, 0);
    }

    std::vector< double > q(steps.size());
    std::vector< double > d(steps.size());
    for (unsigned int i = 0 ; i < steps.size() ; ++i) {
        nvis::vec2 dist = steps[i] - x0;
        q[i] = (double)(i + 1) / dist[0] * metric.width();
        d[i] = metric.distance(x0, steps[i]);
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));

    // bool debug = false;
    // if (drand48()<0.005) {
    //  std::cerr << "period assessment:\n";
    //  std::cerr << "distances = \n";
    //  for (int i=0 ; i<steps.size() ;++i) {
    //      std::cerr << i << ": " << d[i] << ", q = " << q[i] << std::endl;
    //  }
    //  std::cerr << "min distance = " << min << std::endl;
    //  debug = true;
    // }

    double total_q = 0;
    double total_w = 0;
    // if (debug) {
    //  std::cerr << "resulting weights = \n";
    // }
    for (int i = 0 ; i < steps.size() ; ++i) {
        double w = min / fabs(d[i]);
        total_w += w;
        total_q += w * q[i];
        // if (debug) {
        //  std::cerr << i << ": " << w << std::endl;
        // }
    }
    // std::cerr << "returned value = " << total_q << "/" <<  total_w << " = "
    //           << total_q / total_w << std::endl;

    if (total_q / total_w < 0) {
        // std::cerr << "invalid period: " << total_q << " / " << total_w << std::endl;
        // std::cerr << "setting to 0\n";
        // if (x0[1]>40)
        nvis::bbox2 box = bounds(steps);
        if (box.max()[1] > 32) {
            output_broken_chain(steps, metric);
        }
        total_q = 0.;
    }

    return std::make_pair(total_q / total_w, 1.);
}

#ifdef FFTW3
std::pair<double, double>
inline period_x_periodic_fourier(const std::vector<nvis::vec2>& steps,
                                 const nvis::vec2& x0, const metric_type& metric)
{
    if (!steps.size()) {
        return std::make_pair(-1, 0);
    }

    std::vector< double > q(steps.size() + 1);
    q[0] = 0;
    for (int i = 0 ; i < steps.size() ; ++i) {
        nvis::vec2 dp = metric.distance(x0, steps[i]);
        q[i+1] = dp[0];
    }
    return std::make_pair(fft(q), 0.);
}
#endif

void connect_by_distance(std::vector<std::vector<int> >& curves,
                         const std::vector<nvis::vec2>& steps,
                         const spurt::default_metric_type& metric);

void connect_and_check(std::vector<std::vector<int> >& curves,
                       const std::vector<nvis::vec2>& steps,
                       const spurt::default_metric_type& metric);

void connect_by_period(std::vector<std::vector<int> >& curves,
                       const std::vector<nvis::vec2>& steps,
                       const spurt::default_metric_type& metric);

inline double average_distance(const std::vector<nvis::vec2>& steps, unsigned int period,
                               const spurt::default_metric_type& metric)
{
    if (period >= steps.size()) {
        return std::numeric_limits<double>::max();
    }
    double c = 0;
    for (int i = 0 ; i < steps.size() - period ; ++i) {
        c += metric.distance(steps[i], steps[i+period]);
    }
    return c / (double)(steps.size() - period);
}

struct distance_profile {
    distance_profile(const std::vector<nvis::vec2>& steps, unsigned int period,
                     const spurt::default_metric_type& metric) : _p(period) {
        _mean = _median = std::numeric_limits<double>::max();
        _dist.resize(steps.size() - _p);
        double sum = 0;
        for (int i=0 ; i+_p<steps.size() ; ++i) {
            _dist[i] = metric.distance(steps[i], steps[i+period]);
            sum += _dist[i];
        }
        _mean = sum/(double)(_dist.size());
        std::vector<double> tmp(_dist.begin(), _dist.end());
        std::sort(tmp.begin(), tmp.end());
        _median = tmp[tmp.size()/2];
    }

    int _p;
    std::vector<double> _dist;
    double _mean, _median;
};

struct Lt_dist_profile {
    bool operator()(const distance_profile& dp0, const distance_profile& dp1) const {
        if (!dp0._dist.size()) {
            return false;
        } else if (!dp1._dist.size()) {
            return true;
        }

        int N = std::min(dp0._dist.size(), dp1._dist.size());
        int lt_count = 0;
        for (int i=0 ; i<N ; ++i) {
            if (dp0._dist[i] < dp1._dist[i]) {
                ++lt_count;
            }
        }
        if (lt_count > N-lt_count) {
            return true;
        } else if (lt_count < N-lt_count) {
            return false;
        } else {
            return (dp0._median < dp1._median);
        }
    }
};

}

#endif
