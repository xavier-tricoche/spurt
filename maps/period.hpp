#ifndef __XAVIER_PERIOD_HPP__
#define __XAVIER_PERIOD_HPP__

#include <vector>
#include <iostream>
#include <map>
#include <math/fixed_vector.hpp>
#include <poincare/metric.hpp>
#include "cohen_sutherland.hpp"

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

template<typename T1, typename T2>
struct Lt_pair_lexicographic {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1) {
        return (p0.first < p1.first) ||
               (p0.first == p1.first && p0.second < p1.second);
    }
};

template<typename T1, typename T2>
struct Lt_pair_second {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1) {
        return (p0.second < p1.second);
    }
};

}

namespace spurt {

void clip(std::list<std::pair<nvis::vec2, nvis::vec2> >& clipped,
          const nvis::vec2& _x, const nvis::vec2& _y,
          const map_metric& metric)
{
    const map_metric::bounds_type& bounds = metric.bounds();
    
    nvis::vec2 xy = metric.displacement(_x, _y);
    nvis::vec2 x = metric.modulo(_x);
    nvis::vec2 y = metric.modulo(_y);
    clipped.clear();
    if (bounds.inside(x+xy)) {
        clipped.push_back(std::make_pair(x, x+xy));
    } else {
        std::pair<nvis::vec2, nvis::vec2> tmp1, tmp2;
        cohen_sutherland::clip(tmp1, x, x+xy, bounds);
        cohen_sutherland::clip(tmp2, y, y-xy, bounds);
        clipped.push_back(tmp1);
        clipped.push_back(tmp2);
    }
}

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

namespace {
struct Lt_fabs {
    bool operator()(double a, double b) {
        return fabs(a) < fabs(b);
    }
};
}

inline void output_broken_chain(const std::vector<nvis::vec2>& steps,
                                const map_metric& metric)
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
                           const map_metric& metric)
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
dist_based_x_period(const std::vector<nvis::vec2>& steps, const map_metric& metric)
{
    // NB: this method has quadratic complexity...
    
    typedef std::pair<unsigned int, unsigned int>  edge_id_type;
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
    //     std::cerr << "i=" << i << ", j=" << j << ", dx="
    //     << steps[j][0] - steps[i][0] << ", q=" << q << '\n';
    // }
    return q;
}

std::pair<double, double>
inline period_x_periodic(const std::vector<nvis::vec2>& steps, const map_metric& metric)
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
                         const nvis::vec2& x0, const map_metric& metric)
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
    //     std::cerr << "period assessment:\n";
    //     std::cerr << "distances = \n";
    //     for (int i=0 ; i<steps.size() ;++i) {
    //         std::cerr << i << ": " << d[i] << ", q = " << q[i] << std::endl;
    //     }
    //     std::cerr << "min distance = " << min << std::endl;
    //     debug = true;
    // }
    
    double total_q = 0;
    double total_w = 0;
    // if (debug) {
    //     std::cerr << "resulting weights = \n";
    // }
    for (int i = 0 ; i < steps.size() ; ++i) {
        double w = min / fabs(d[i]);
        total_w += w;
        total_q += w * q[i];
        // if (debug) {
        //    std::cerr << i << ": " << w << std::endl;
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
                                 const nvis::vec2& x0, const map_metric& metric)
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


}

#endif













































