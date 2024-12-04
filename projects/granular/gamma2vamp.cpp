#include <iostream>
#include <math.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "USAGE: " << argv[0] << " <gamma>\n";
        exit(-1);
    }
    double gamma = atof(argv[1]);
    double omega = 2.*M_PI*7.5;
    double amp = 9.81*gamma/(omega*omega);
    double vamp = omega*amp;
    std::cout << amp << '\t' << vamp << '\n';
    return 0;
}