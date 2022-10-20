#include "XY2LatLong.hpp"

int main(int argc, char* argv[]) {
    spurt::XY2LatLong converter(argv[1]);
    
    while(true) {
        std::cout << "Enter Lat/Long: " << '\n';
        double la, lo;
        std::cin >> la >> lo;
        std::cout << "corresponding coordinates are " << converter.toXY(nvis::vec2(la, lo)) << '\n';
    }
    exit(0);
}