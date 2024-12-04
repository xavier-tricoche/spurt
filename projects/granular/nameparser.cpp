#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, const char* argv[])
{
    if (argc != 2) {
        std::cerr << "USAGE: " << argv[0] << " <filename>\n";
        return -1;
    }
    
    std::string name(argv[1]);
    int n = name.size();
    
    size_t ext = name.rfind(".nrrd");
    if (ext == std::string::npos) {
        std::cerr << "input file name does not contain .nrrd extension\n";
        return -1;
    }
    
    size_t tpos = name.rfind("-t=");
    if (tpos == std::string::npos) {
        std::cerr << "input file name does not contain time information\n";
        return -1;
    } else if (tpos + 3 >= ext) {
        std::cerr << "input file name has wrong format\n";
        return -1;
    }
    
    size_t tlength = ext - tpos - 3;
    char tchar[tlength];
    name.copy(tchar, tlength, tpos+3);
    float t = atof(tchar);
    std::cout << "time = " << t << '\n';
    
    return 0;
}



