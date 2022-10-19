#include <iostream>
#include <teem/limn.h>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cout << "USAGE: " << argv[0] << " <input.lmpd> <output.vtk" << std::endl;
        return -1;
    }
    
    limnPolyData* lpd = limnPolyDataNew();
    limnPolyDataReadLMPD(lpd, fopen(argv[1], "r"));
    
    FILE* file = fopen(argv[2], "w");
    limnPolyDataWriteVTK(file, lpd);
    
    return 0;
}

