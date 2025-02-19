#include <vector>
#include <set>
#include <algorithm>

namespace spurt {

template<typename Iter_, typename Val_ = typename Iter_::value_type>
Val_ min(Iter_ begin, Iter_ end) {
	return *std::min_element(begin, end);
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
Val_ max(Iter_ begin, Iter_ end) {
	return *std::max_element(begin, end);
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
size_t minrank(Iter_ begin, Iter_ end) {
	return std::distance(begin, std::min_element(begin, end));
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
size_t maxrank(Iter_ begin, Iter_ end) {
	return std::distance(begin, std::max_element(begin, end));
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
std::pair<Val_, Val_> minmax(Iter_ begin, Iter_ end) {
	auto mm = std::minmax_element(begin, end);
	return std::make_pair(*mm.first, *mm.second);
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
std::pair<double, double> meanvariance(Iter_ begin, Iter_ end)
{
    if (begin == end) return std::pair<double, double>(0, 0);
	size_t len = std::distance(begin, end);
    double mean = 0;
    for (typename std::vector<Val_>::const_iterator it = begin ; it != end ; ++it) {
        mean += *it;
    }
    mean /= (double)len;
    
    double variance = 0;
    for (typename std::vector<Val_>::const_iterator it = begin ; it != end ; ++it) {
        double d = (double)(*it) - mean;
        variance += d * d;
    }
    variance /= (double)(len-1);
    
    return std::make_pair(mean, variance);
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
Val_ median(Iter_ begin, Iter_ end)
{
	std::vector<Val_> copy(begin, end);
	if (copy.empty()) {
		std::cerr << "WARNING: attempting to compute median of empty set\n";
		return Val_(0);
	}
	else if (copy.size() == 1) return copy[0];

	std::sort(copy.begin(), copy.end());
	if (copy.size() % 2) {
		// n=3 -> size/2 = 1
		return copy[copy.size()/2];
	}
	else {
		// n=4 -> size/2-1 = 1, size/2 = 2 -> { 0, 1, 2, 3 } -> 1, 2
		Val_ low = copy[copy.size()/2-1];
		Val_ high = copy[copy.size()/2];
		return 0.5*(low+high);
	}
}

template<typename Iter_, typename Val_ = typename Iter_::value_type>
std::vector<double> percentiles(Iter_ begin, Iter_ end, int n=11)
{
	std::vector<double> copy(begin, end);
	if (copy.size() <= 1) return copy; 

	std::sort(copy.begin(), copy.end());
	if (copy.size() <= n) return copy;
	auto stride = std::div(copy.size(), n-1);
	std::vector<size_t> steps;
	int q = stride.quot;
	int r = stride.rem;
	std::cout << "std::div(" << copy.size() << ", " << n-1 << ") = " << q << ", " << r << '\n';
	steps.push_back(0);
	for (int i=1; i<n-1; ++i, --r) {
		int s = steps.back() + q;
		std::cout << "i=" << i;
		int olds = s;
		if (r > 0) {
			++s;
		}
		std::cout << ", s=" << s << " (r=" << r << ", " << olds << ")\n";
		steps.push_back(s);
	}
	steps.push_back(copy.size()-1);
	std::vector<double> res;
	std::for_each(steps.begin(), steps.end(), [&](int s) {
		res.push_back(copy[s]);
	});
	return res;
}

template<typename T> 
T max(const std::vector<T>& values) {
	return spurt::max(values.begin(), values.end());
}

template<typename T> 
T min(const std::vector<T>& values) {
	return spurt::min(values.begin(), values.end());
}

template<typename T>
size_t maxrank(const std::vector<T>& values) {
	return maxrank(values.begin(), values.end());
}

template<typename T>
size_t minrank(const std::vector<T>& values) {
	return minrank(values.begin(), values.end());
}

template<typename T>
std::pair<T, T> minmax(const std::vector<T>& values) {
	return minmax(values.begin(), values.end());
}

template<typename T>
std::pair<double, double> meanvariance(const std::vector<T>& values)
{
    if (values.empty()) return std::pair<double, double>(0, 0);
    double mean = 0;
    for (typename std::vector<T>::const_iterator it = values.begin() ; it != values.end() ; ++it) {
        mean += *it;
    }
    mean /= (double)values.size();
    
    double variance = 0;
    for (typename std::vector<T>::const_iterator it = values.begin() ; it != values.end() ; ++it) {
        double d = (double)(*it) - mean;
        variance += d * d;
    }
    variance /= (double)(values.size()-1);
    
    return std::make_pair(mean, variance);
}

template<typename T>
T median(const std::vector<T>& values) 
{
	return median(values.begin(), values.end());
}

template<typename T>
std::vector<double> percentiles(const std::vector<T>& values, int n) 
{
	return percentiles(values.begin(), values.end(), n);
}

template<typename T>
struct histogram_info {
	std::vector<size_t> bins;
	size_t nbins;
	T min_val;
	T max_val;
	double delta;
};

template<typename Iter, typename T = typename Iter::value_type>
histogram_info<T> histo(Iter begin, Iter end, size_t nbins,
					    T min_ = std::numeric_limits<T>::min(), 
					    T max_ = std::numeric_limits<T>::max())
{
	auto mM = std::minmax_element(begin, end);
	T actual_min = *mM.first;
	T actual_max = *mM.second;
	if (min_ == std::numeric_limits<T>::min()) {
		min_ = actual_min;
	}
	else {
		min_ = std::max(min_, actual_min);
	}
	if (max_ == std::numeric_limits<T>::max()) {
		max_ = actual_max;
	}
	else {
		max_ = std::min(max_, actual_max);
	}

	double delta = (max_ - min_) / static_cast<double>(nbins);
	histogram_info<T> r;
	r.min_val = min_;
	r.max_val = max_;
	r.delta = delta;
	r.nbins = nbins;
	r.bins.resize(nbins, 0);

	std::for_each(begin, end, [&](T v) {
		double i = (v-min_) / delta;
		r.bins[std::floor(i)]++;
	});

	return r;
}

template<typename Iter, typename T = typename Iter::value_type>
double mode(Iter begin, Iter end, size_t nbins=1024,
			T min_ = std::numeric_limits<T>::min(),
			T max_ = std::numeric_limits<T>::max())
{
	histogram_info<T> h = histo(begin, end, nbins, min_, max_);
	size_t maxbin = maxrank(h.bins);

	return h.min_val + maxbin*h.delta + 0.5*h.delta; 
}

template<typename T>
double mode(const std::vector<T>& values, size_t nbins=1024,
			T min_ = std::numeric_limits<T>::min(),
			T max_ = std::numeric_limits<T>::max())
{
	return mode(values.begin(), values.end(), nbins, min_, max_);
}

template<typename T>
histogram_info<T> histo(const std::vector<T>& values, size_t nbins,
					    T min_ = std::numeric_limits<T>::min(), 
					    T max_ = std::numeric_limits<T>::max())
{
	return histo(values.begin(), values.end(), min_, max_);
}

}



