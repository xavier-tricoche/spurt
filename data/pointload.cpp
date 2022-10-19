#include <iostream>
#include <teem/nrrd.h>
#include <limits>
#include <math/fixed_vector.hpp>

const double load = 1;
const double poisson_ratio = 0.3;

int main(int argc, char* argv[])
{
    if (argc != 8) {
        std::cout << "USAGE: " << argv[0] << " <length> <width> <depth> <Nx> <Ny> <Nz> <output>"
                  << std::endl;
        return -1;
    }
    
    float length = atof(argv[1]);
    float width  = atof(argv[2]);
    float depth  = atof(argv[3]);
    vec3 origin(-0.5*length, -0.5*width, -depth);
    
    unsigned int Nx = atoi(argv[4]);
    unsigned int Ny = atoi(argv[5]);
    unsigned int Nz = atoi(argv[6]);
    vec3 spacing(length / (Nx - 1), width / (Ny - 1), depth / (Nz - 1));
    
    int i, j, k;
    double tensor[9];
    unsigned int numPts;
    double P, twoPi, rho, rho2, rho3, rho5, nu;
    double x, x2, y, y2, z, z2, rhoPlusz2, zPlus2rho, txy, txz, tyz;
    vec3 s, center[2];
    
    std::cout << "Computing point load stress tensors" << std::endl;
    
    //
    // Initialize self; create output objects
    //
    numPts = Nx * Ny * Nz;
    
    double* data = (double*)calloc(Nx * Ny * Nz * 7, sizeof(double));
    
    //
    // Compute the location of the load
    //
    center[0] = origin + vec3(0.25 * length, 0.5 * width, depth);
    center[1] = origin + vec3(0.75 * length, 0.5 * width, depth);
    
    for (unsigned int n = 0 ; n < 2 ; ++n) {
        // Traverse all points evaluating implicit function at each point. Note that
        // points are evaluated in local coordinate system of applied force.
        //
        twoPi = 2.0 * M_PI;
        P = -load;
        int pointCount = 0;
        for (k = 0; k < Nz; k++) {
            z = center[n][2] - (origin[2] + k * spacing[2]);
            for (j = 0; j < Ny; j++) {
                y = center[n][1] - (origin[1] + j * spacing[1]);
                for (i = 0; i < Nx; i++) {
                    x = (origin[0] + i * spacing[0]) - center[n][0];
                    rho = sqrt(x * x + y * y + z * z);//in local coordinates
                    double* tensor = &data[7 * (i + Nx * (j + Ny * k))];
                    tensor[0] = 1;
                    if (rho < 1.0e-10) {
                        std::cerr << "Attempting to set singularity, resetting" << std::endl;
                        
                        tensor[1] += std::numeric_limits<double>::max(); // Component(0,0)
                        tensor[4] += std::numeric_limits<double>::max(); // Component(1,1);
                        tensor[6] += std::numeric_limits<double>::max(); // Component(2,2);
                        
                        tensor[2] += 0.;
                        tensor[3] += 0.;
                        tensor[5] += 0.;
                        
                        pointCount++;
                        continue;
                    }
                    
                    rho2 = rho * rho;
                    rho3 = rho2 * rho;
                    rho5 = rho2 * rho3;
                    nu = (1.0 - 2.0 * poisson_ratio);
                    x2 = x * x;
                    y2 = y * y;
                    z2 = z * z;
                    rhoPlusz2 = (rho + z) * (rho + z);
                    zPlus2rho = (2.0 * rho + z);
                    
                    // normal stresses
                    s[0] = P / (twoPi * rho2) * (3.0 * z * x2 / rho3 - nu * (z / rho - rho / (rho + z) +
                                                 x2 * (zPlus2rho) / (rho * rhoPlusz2)));
                    s[1] = P / (twoPi * rho2) * (3.0 * z * y2 / rho3 - nu * (z / rho - rho / (rho + z) +
                                                 y2 * (zPlus2rho) / (rho * rhoPlusz2)));
                    s[2] = 3.0 * P * z2 * z / (twoPi * rho5);
                    
                    //shear stresses - negative signs are coordinate transformations
                    //that is, equations (in text) are in different coordinate system
                    //than volume is in.
                    txy = -(P / (twoPi * rho2) * (3.0 * x * y * z / rho3 -
                                                  nu * x * y * (zPlus2rho) / (rho * rhoPlusz2)));
                    txz = -(3.0 * P * x * z2 / (twoPi * rho5));
                    tyz = 3.0 * P * y * z2 / (twoPi * rho5);
                    
                    tensor[1] += s[0];  // Component(0,0);
                    tensor[4] += s[1];  // Component(1,1);
                    tensor[6] += s[2];  // Component(2,2);
                    tensor[2] += txy;   // Component(0,1);  real symmetric matrix
                    tensor[3] += txz;   // Component(0,2);
                    tensor[5] += tyz;   // Component(1,2);
                    
                    pointCount++;
                }
            }
        }
    }
    
    // export in NRRD file
    char* err;
    Nrrd* nout = nrrdNew();
    if (nrrdWrap_va(nout, data, nrrdTypeDouble, 4, 7, Nx, Ny, Nz)) {
        err = biffGetDone(NRRD);
        std::cerr << "trouble wrapping data volume: " << err << std::endl;
        free(err);
        return -1;
    }
    
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoLabel, "Conf;Txx;Txy;Txz;Tyy;Tyz;Tzz", "x", "y", "z");
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterUnknown, nrrdCenterNode, nrrdCenterNode, nrrdCenterNode);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, AIR_NAN, spacing[0], spacing[1], spacing[2]);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoMin, AIR_NAN, origin[0], origin[1], origin[2]);
    
    if (nrrdSave(argv[7], nout, NULL)) {
        err = biffGetDone(NRRD);
        std::cerr << "trouble saving nrrd file to disk: " << err << std::endl;
        free(err);
        return -1;
    }
    
    return 0;
}























