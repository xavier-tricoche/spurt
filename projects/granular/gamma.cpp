#include <iostream>
#include <math.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    // USAGE: argv[0] <gamma>
    if (argc != 3) {
        std::cerr << "USAGE: " << argv[0] << " <gamma> <freq>\n";
        exit(-1);
    }
    double gamma = atof(argv[1]);
    double f = atof(argv[2]);
    
    // Γ ≡ aω2/g, ω = 2πf
    double omega = 2.*M_PI*f;
    const double g = 9.81;
    double vamp = gamma*g/omega;
    double a = gamma*g/(omega*omega);
    std::cout << "Gamma=" << gamma << ", frequency=" << f
              << ", omega=" << omega << ", amp=" << a
              << ", vamp=" << vamp << '\n';
              
    return 1;
}