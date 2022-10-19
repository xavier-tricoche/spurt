#include <fstream>
#include <iostream>


int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "USAGE: " << argv[0] << " <in.vtk> <out.obj>\n";
        exit(-1);
    }
    
    std::fstream in(argv[1], std::ios::in);
    if (!in) {
        std::cerr << "ERROR: " << argv[1] << ": no such file\n";
        exit(-1);
    }
}
