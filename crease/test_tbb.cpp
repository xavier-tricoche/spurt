#include <map>
#include <thread>
#include <iostream>

// TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>

int main(int argc, char** argv) {
    std::map<float, float> amap;
    typedef std::map<float, float>::value_type value_type;

    for (long i=0; i<1000000; ++i) {
        auto k = 1./float(i+1);
        auto v = exp(k);
        amap[k] = v;
    }
    std::cout << "map initialized\n";

    std::vector<float> keys;
    for (auto kv : amap) {
        keys.push_back(kv.first);
    }

    std::cout << "keys copied\n";

    std::atomic<size_t> counter = 0;
    std::atomic<float> keysum = 0;
    std::atomic<float> valsum = 0;
    tbb::parallel_for(tbb::blocked_range<long>(0, keys.size()),
                      [&](tbb::blocked_range<long> r)
    {
        for (auto n=r.begin(); n!=r.end(); ++n)
        {
            ++counter;
            auto k = keys[n];
            valsum = valsum + amap[k];
            keysum = keysum + k;
        }
    });

    std::cout << "there were " << counter << " iterations. Sum of keys is " << keysum << ", sum of values is " << valsum << '\n';

    return 0;
}