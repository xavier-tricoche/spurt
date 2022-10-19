#include <vector>
#include <set>

namespace xavier
{

template<typename T>
std::pair<double, double> meanvariance(const std::vector<T>& values)
{
	if (values.empty()) {
		return std::pair<double, double>(0, 0);
	}
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
	variance /= (double)values.size();

	return std::make_pair(mean, variance);
}

}



