#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <blaze/Blaze.h>

template<typename T>
T r() {
    return static_cast<T>(std::rand())/static_cast<T>(RAND_MAX);
}

typedef blaze::StaticVector<double, 2UL> vec2;
typedef blaze::StaticVector<double, 3UL> vec3;
typedef blaze::StaticVector<double, 4UL> vec4;
typedef blaze::StaticVector<float, 3UL> fvec3;
typedef blaze::StaticVector<int, 3UL> ivec3;

int main(int argc, char* argv[]) {

    vec3 d[4];
    fvec3 f[4];
    ivec3 j[4];

    std::srand(time(nullptr));
    for (int i=0; i<4; ++i) {
        double _a = r<double>();
        double _b = r<double>();
        double _c = r<double>();
        d[i] = vec3({_a, _b, _c});
        f[i] = fvec3({r<float>(), r<float>(), r<float>()});
        j[i] = ivec3({int(r<float>()*100), int(r<float>()*100), int(r<float>()*100)});
    }
    vec2 a({r<double>(), r<double>()});
    vec4 b({r<double>(), r<double>(), r<double>(), r<double>()});

    vec3 l = d[0] + d[1];
    std::cout << "l=" << l << '\n';
    vec3 s = d[2] + d[3];
    if (blaze::any(s>0.5)) std::cout << "any TRUE\n";
    else std::cout << "any is FALSE\n";
    fvec3 g = f[0]*f[1];
    fvec3 h = f[2]/f[3];
    std::cout << "the norm of " << f[0] << " * " << f[1] <<  " (=" << g << ") is " << norm(f[0]*f[1]) << '\n';
    std::cout << "the norm of " << f[2] << " / " << f[3] << " (=" << h << ") is " << norm(f[2]/f[3]) << '\n';

    vec3 k(f[0]);
    std::cout << "Initializing a double vector with a float vector produced " << k << '\n';

    std::cout << "The sum of double vector " << d[0] << " and a float vector " << f[0] << " is " << d[0]+f[0] << '\n';
    std::cout << "The sum of a float vector " << f[1] << " and an int vector " << j[0] << " is " << f[1]+j[0] << '\n';
    std::cout << "The product of a double vector " << d[1] << " and an int vector " << j[1] << " is " << d[1] * j[1] << '\n';

    std::sort(d, d+4);
    std::cout << "after sorting, the vectors are:\n";
    for (int i=0; i<4; ++i) std::cout << d[i] << '\n';

    d[2][1] *= -1;
    std::cout << "the absolute value of " << d[2] << " is " << abs(d[2]) << '\n';

    std::cout << "the square of " << d[3] << " is " << square(d[3]) << '\n';


    return 0;
}