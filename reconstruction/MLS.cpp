#ifdef _OPENMP
#include <omp.h>
#endif
#include <mls.h>
#include <kdtree++/kdtree.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <teem/nrrd.h>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <boost\math\constants\constants.hpp>
#include <math.h>

// parameters
size_t             npts;
unsigned int     resolution;
char*             basename;
double             radius;

void initialize(int argc, char* argv[])
{
    hestOpt *hopt = NULL;
    hestParm *hparm;
    airArray *mop;
    char *me;

    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "n",    "# points",        airTypeSize_t,     0,     1,    &npts,            "10000",        "number of sample points");
    hestOptAdd(&hopt, "s",    "size",            airTypeInt,     0,     1,    &resolution,    "512",            "size of (square) output image");
    hestOptAdd(&hopt, "r",    "radius",        airTypeDouble,     0,     1,    &radius,        "0.05",            "radius of MLS fit");
    hestOptAdd(&hopt, "o",    "output",        airTypeString,     1,     1,    &basename,        NULL,            "output base name");

    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Reconstruct Franke test function using MLS",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T, int N>
class point {
public:
    typedef T                             value_type;
    typedef nvis::fixed_vector<T, N>    vector_type;
    typedef size_t                        index_type;

    point() : __i(), __v() {}    
    point(const index_type& i, const vector_type& v) : __i(i), __v(v) {}
    point(const point& p) : __i(p.__i), __v(p.__v) {}

    const vector_type& pos() const {
        return __v;
    }

    vector_type& pos() {
        return __v;
    }

    index_type idx() const {
        return __i;
    }

    index_type& idx() {
        return __i;
    }

    value_type distance_to(const point& p) const {
        return norm(__v - p.__v);
    }

    value_type operator[](size_t i) const {
        return __v[i];
    }

private:
    index_type    __i;
    vector_type __v;
};


const double pi = M_PI;

inline double sqr(double x){
    return x*x;
}

inline double cube(double x){
    return x*x*x;
}




//makes a random double
static double drand(){
    //return rand() / (RAND_MAX + 1.0);
    return ((double)(rand() % 100000 - 50000.0))/50000.0;
}

int main(int argc, char* argv[]) {
    const int dim = 3;
    const int nrhs = 1;


    typedef point<double, dim>                  point_type;
    typedef nvis::fixed_vector<double, nrhs>      value_type;
    typedef point_type::vector_type                vector_type;
    
    std::vector<value_type>                     values;
    initialize(argc, argv);
    
    std::cerr << "radius = " << radius << std::endl;

    KDTree::KDTree<dim, point_type> tree;

    // //finds points in the franke function
    // int npts = /*atoi(argv[1])*/ 10000;
    double *sdata = (double*)calloc(npts*3, sizeof(double));

    values.resize(npts);
    for(int i = 0; i < npts; i++){
        double x = drand();
        double y = drand();
        double z = drand();
        double x1 = MarschnerLobb(x,y,z);
        
        vector_type pnt(x,y,z);
        values[i][0] = x1;
        tree.insert(point_type(i, pnt));
        sdata[4*i+0] = x;
        sdata[4*i+1] = y;
        sdata[4*i+2] = z;
        sdata[4*i+3] = x1;
    }
    std::cout << "Making -samples.nrrd" << std::endl;
    std::vector<size_t> ssize(2);
    ssize[0] = 4;
    ssize[1] = npts;
    std::string sample_name = basename;
    sample_name += "-samples.nrrd";
    spurt::writeNrrd(sdata, sample_name, nrrdTypeDouble, ssize);
    std::cout << "Made -samples.nrrd" << std::endl;

    int n_runs = resolution*resolution*resolution;
    double *data = (double*)calloc(n_runs, sizeof(double));

    MLS::polynomial_fit<value_type, vector_type, dim> pfit = MLS::polynomial_fit<value_type, vector_type, dim>(dim, 0, nrhs);
    clock_t init, final;
    double rsearchtime=0;
    init=clock();
    int count = 0;
    float *ground_truth = (float*)calloc(resolution*resolution*resolution, sizeof(float));
    std::cout << "Starting Production" << std::endl;
    for(int n = (-1*((int)resolution))/2; n < (int)resolution/2 ; n++){
        for(int m = (-1*((int)resolution))/2; m < (int)resolution/2 ; m++){
            for(int l = (-1*((int)resolution))/2; l < (int)resolution/2 ; l++){
                vector_type p;
                p[0] = (double)n/((double)resolution/2);
                p[1] = (double)m/((double)resolution/2);
                p[2] = (double)l/((double)resolution/2);
                


                ground_truth[count] = MarschnerLobb(p[0], p[1], p[2]);
                
                std::cout << p[0] << ", " << p[1] << ", " << p[2] << "    " << MarschnerLobb(p[0], p[1], p[2]) << std::endl;
            
                point_type pt(0, p);

                // in radius box find
                std::vector<point_type > in_cube;
                clock_t sinit, sfinal;
                sinit=clock();
                //tree.find_within_range(pt, radius, std::back_inserter(in_cube)); 
                sfinal=clock()-sinit; 
                rsearchtime += (double)sfinal / ((double)CLOCKS_PER_SEC);
                // eliminate points not actually within radius
                std::vector<point_type> valid;
                for(int i = 0; i < in_cube.size(); i++) {
                    if (nvis::norm(p-in_cube[i].pos()) < radius) {
                        valid.push_back(in_cube[i]);
                    }
                }
                Eigen::MatrixXd fitted_coef;
                std::vector<value_type> rvalues;

                std::vector<vector_type> rpoints;
                if(valid.empty()){
                    //std::cerr << "No Points within Radius" << std::endl;
                    continue;
                }
                // std::cerr << in_cube.size() << " / " << valid.size() << " found\n";

                //if(m*n*l <= 1500000){
                //    std::cout << "almost there!" << std::endl;
                //}

                int c_size = valid.size();
                rpoints.resize(c_size);
                rvalues.resize(c_size);

                for(int i = 0; i < c_size; i++){
                    rpoints[i] = valid[i].pos();
                    rvalues[i] = values[valid[i].idx()];
                }
                //int prec = pfit(fitted_coef, rpoints, rvalues, p, radius);  
                //int a = fitted_coef.size();
                data[count] = 0/*fitted_coef(0)*/;
                count ++;
                //std::cout << count << std::endl;
                //std::cout << count <</* ", " << fitted_coef(0) << ", " << franke(p[0],p[1]) <<" (" << p[0] << ", " << p[1] << ") " << in_cube.size() << */std::endl;
            }
        }
    }

    //for(unsigned int n=0 ; n < resolution ; n++){
    //    for(unsigned int m = 0 ; m < resolution ; m++){
    //        vector_type p;
    //        p[0] = (double)n/(double)resolution;
    //        p[1] = (double)m/(double)resolution;
    //        
    //        ground_truth[m*resolution+n] = /*f7*//*f10*/franke(p[0], p[1]);
    //        
    //        point_type pt(0, p);

    //        // in radius box find
    //        std::vector<point_type > in_cube;
    //        clock_t sinit, sfinal;
    //        sinit=clock();
    //        tree.find_within_range(pt, radius, std::back_inserter(in_cube)); 
    //        sfinal=clock()-sinit; 
    //        rsearchtime += (double)sfinal / ((double)CLOCKS_PER_SEC);
    //        // eliminate points not actually within radius
    //        std::vector<point_type> valid;
    //        for(int i = 0; i < in_cube.size(); i++) {
    //            if (nvis::norm(p-in_cube[i].pos()) < radius) {
    //                valid.push_back(in_cube[i]);
    //            }
    //        }
    //        Eigen::MatrixXd fitted_coef;
    //        std::vector<value_type> rvalues;

    //        std::vector<vector_type> rpoints;
    //        if(valid.empty()){
    //            // std::cerr << "No Points within Radius" << std::endl;
    //            continue;
    //        }
    //        // std::cerr << in_cube.size() << " / " << valid.size() << " found\n";

    //        int c_size = valid.size();
    //        rpoints.resize(c_size);
    //        rvalues.resize(c_size);

    //        for(int i = 0; i < c_size; i++){
    //            rpoints[i] = valid[i].pos();
    //            rvalues[i] = values[valid[i].idx()];
    //        }
    //        int prec = pfit(fitted_coef, rpoints, rvalues, p, radius);  
    //        int a = fitted_coef.size();

    //        data[m*resolution+n] = fitted_coef(0);
    //        //std::cout << count << ", " << fitted_coef(0) << ", " << franke(p[0],p[1]) <<" (" << p[0] << ", " << p[1] << ") " << in_cube.size() << std::endl;
    //        count ++;
    //    }
    //}
    final=clock()-init;
    std::cout << "Algorithm time : " << (double)final / ((double)CLOCKS_PER_SEC) << " for " << resolution*resolution << " runs, on " << npts << " points"  <<std::endl;
    std::cout << rsearchtime << " seconds spent finding points" << std::endl;
    std::cout << pfit.svdtime << " seconds spent doing SVD" << std::endl;


    // NRRD file storage
    std::vector<size_t> size(dim);
    std::vector<double> spacing(dim);
    spacing[0] = 1./(double)resolution;
    spacing[1] = 1./(double)resolution;
    spacing[2] = 1./(double)resolution;
    size[0] = resolution;
    size[1] = resolution;
    size[2] = resolution;
    std::string output_name = basename;
    output_name += "-results.nrrd";
    spurt::writeNrrd(data, output_name, nrrdTypeDouble, size, spacing);
    std::cerr << "output results.NRRD file exported\n";

    output_name = basename;
    output_name += "-ml.nrrd";
    spurt::writeNrrd(ground_truth, output_name, nrrdTypeFloat, size, spacing);
    std::cerr << "output ml.NRRD file exported\n";
    return 0;
}
