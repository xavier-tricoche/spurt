#ifndef __XAVIER_READ_SPARSE__
#define __XAVIER_READ_SPARSE__

#include <image/nrrd_wrappper.hpp>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Triplet;


int main(int argc, char* argv[]) {
    
    std::string name_ids, name_coef;
    name_ids = argv[1];
    name_coef = argv[2];
    
    Nrrd* nin_ids = xavier::nrrd_utils::readNrrd(name_ids);
    Nrrd* nin_coef = xavier::nrrd_utils::readNrrd(name_ids);
    
    int *ids = static_cast<int *>(nin_ids->data);
    size_t n_ids = nin_ids->axis[0].size;

    std::vector<Triplet> triplets;
    size_t row=0;
    for (size_t i=0 ; i<n_ids ; ++i) {
        unsigned int n_cols = ids[i];
        for (unsigned int j=0; j<n_cols; ++j) {
            size_t col=ids[i+j];
            triplets.push_back(Triplet(row, col, 0.));
        }
        i += n_cols;
    }
    
    std::cout << "there are " << triplets.size() << " non-zero terms\n";
    
    nrrdNuke(nin_ids);
    nrrdNuke(nin_coefs);
    
    return 0;
}


#endif // __XAVIER_READ_SPARSE__
