#include <teem/limn.h>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cout << "USAGE: " << argv[0] << " <input.lmpd> <output.vcl>" << std::endl;
        return -1;
    }
    
    limnPolyData* lpd = limnPolyDataNew();
    limnPolyDataReadLMPD(lpd, fopen(argv[1], "r"));
    
    const float* pos = lpd->xyzw;
    size_t npos = lpd->xyzwNum;
    
    const unsigned int* ids = lpd->indx;
    size_t nids = lpd->indxNum;
    
    const unsigned char* types = lpd->type;
    size_t nprim = lpd->primNum;
    const unsigned int* nvert = lpd->icnt;
    
    std::fstream output(argv[2], std::ios::out);
    unsigned int id = 0;
    for (unsigned int i = 0 ; i < nprim ; id += nvert[i], ++i) {
        if (types[i] == limnPrimitiveLineStrip || types[i] == limnPrimitiveLines) {
            for (unsigned int j = 0 ; j < nvert[i] ; ++j) {
                unsigned int pid = ids[id+j];
                const float* x = &pos[4*pid];
                output << "p " << x[0] << " " << x[1] << " " << x[2] << '\n';
            }
            output << "n\n";
        }
    }
    
    output.close();
}



