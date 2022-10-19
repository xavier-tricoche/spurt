#include <vector>
#include <set>
#include <algorithm>

namespace xavier {
	
template<typename T> 
double max(const std::vector<T>& values) {
	if (values.empty()) return 0;
	else return *std::max_element(values.begin(), values.end());
}

template<typename T> 
double min(const std::vector<T>& values) {
	if (values.empty()) return 0;
	else return *std::min_element(values.begin(), values.end());
}

template<typename T>
size_t maxrank(const std::vector<T>& values) {
	if (values.empty()) return 0;
	else return std::distance(values.begin(), std::max_element(values.begin(), values.end()));
}

template<typename T>
size_t minrank(const std::vector<T>& values) {
	if (values.empty()) return 0;
	else return std::distance(values.begin(), std::min_element(values.begin(), values.end()));
}

template<typename T>
std::pair<double, double> minmax(const std::vector<T>& values) {
	if (values.empty()) return std::pair<double, double>(0,0);
	else {
		std::pair<double, double> r;
		auto mm = std::minmax_element(values.begin(), values.end());
		return std::pair<double, double>(*mm.first, *mm.second);
	}
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

}



