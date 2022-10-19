#include <teem/nrrd.h>
#include <ddsbase.h>
#include <image/nrrd_wrapper.hpp>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "USAGE: " << argv[0] << " <in.pvm> <out.nrrd>\n";
        return -1;
    }
    
    unsigned int width, height, depth, components;
    float scalex, scaley, scalez;
    unsigned char* data = readPVMvolume(argv[1], &width, &height, &depth, &components,
                                        &scalex, &scaley, &scalez);
                                        
    std::cerr << "read image " << argv[1] << " which has size "
              << width << " x " << height << " x " << depth
              << " and " << components
              << " char components per pixel\n";
    std::cerr << "scales = " << scalex << " x " << scaley << " x " << scalez << '\n';
    
    
    Nrrd* nout = nrrdNew();
    
    if (nrrdWrap_va(nout, data, nrrdTypeShort, 3, width, height, depth)) {
        std::cerr << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, scalex, scaley, scalez);
    if (nrrdSave(argv[2], nout, NULL)) {
        std::cerr << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    
    return 0;
}




