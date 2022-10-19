#include <kdtree++/kdtree.hpp>
#include <kdtree/kdtree.hpp>
#include <iostream>
#include <util/timer.hpp>
#include <Eigen/Core>
#include <data/indexed_point.hpp>
#include <sfcnn.hpp>

static size_t sizes[] = { 100, 1000, 10000, 100000, 1000000, 10000000 };

template<int dim>
void compare()
{
    typedef Eigen::Matrix<float, dim, 1>        vec_type;
    typedef xavier::indexed_point<float, dim>   point_type;
    typedef kdtree<vec_type, int, dim>          tree_type_1;
    typedef KDTree::KDTree<dim, point_type>     tree_type_2;
    typedef sfcnn<vec_type, dim, float>         tree_type_3;
    
    for (int m=0 ; m<6 ; ++m) {
        int npts = sizes[m];
        
        std::cout << "in " << dim << "D with " << npts << " points:\n";
        
        tree_type_1 tree1;
        tree_type_2 tree2;
        vec_type* all_pts = new vec_type[npts];
        srand48(time(0));
        
        for (int i=0 ; i<npts ; ++i) {
            vec_type    x1;
            point_type  x2;
            for (int d=0 ; d<dim ; ++d) {
                float y = drand48();
                x1[d] = y;
                x2.pos()[d] = y;
            }
            x2.index() = i;
            tree1.add(x1, i);
            tree2.insert(x2);
            all_pts[i] = x1;
        }
        
        nvis::timer init_t;
        tree_type_3 tree3(all_pts, npts);
        std::cerr << "STANN initialization took " << init_t.elapsed() << " s.\n";
        init_t.restart();
        tree1.sort();
        std::cerr << "pointerless kdtree initialization took " << init_t.elapsed() << " s.\n";
        
        double t1 = 0;
        double t2 = 0;
        double t3 = 0;
        double dist1 = 0;
        double dist2 = 0;
        double dist3 = 0;
        nvis::timer t;
        for (int i=0 ; i<100 ; ++i) {
            vec_type    x1;
            point_type  x2;
            for (int d=0 ; d<dim ; ++d) {
                float y = drand48();
                x1[d] = y;
                x2.pos()[d] = y;
            }
            
            std::vector<typename tree_type_1::const_iterator > ans1;
            std::vector<long unsigned int> ans3;
            
            t.restart();
            tree1.find_n_nearest(x1, 20, std::back_inserter(ans1));
            t1 += t.elapsed();
            
            t.restart();
            tree3.ksearch(x1, 20, ans3);
            t3 += t.elapsed();
            
            for (int n=0 ; n<ans1.size() ; ++n) {
                dist1 += (x1 - ans1[n]->first).norm();
            }
            if (ans1.size() < 20) {
                // std::cerr << "kdtree found only " << ans1.size() << " neighbors!\n";
                dist1 += (20-ans1.size()) * sqrt(2);
            }
            
            for (int n=0 ; n<ans3.size() ; ++n) {
                dist3 += (x1 - all_pts[ans3[n]]).norm();
            }
            if (ans3.size() < 20) {
                // std::cerr << "STANN found only " << ans3.size() << " neighbors!\n";
                dist3 += (20-ans3.size()) * sqrt(2);
            }
        }
        
        std::cout  << "\tkdtree: time =" << t1/100. << " s. / dist = " << dist1/2000. << "\n"
                   << "\tSTANN:  time =" << t3/100. << " s. / dist = " << dist3/2000. << "\n\n";
                   
        delete[] all_pts;
    }
}

int main(int argc, char* argv[])
{

    compare<1>();
    compare<2>();
    compare<3>();
    compare<4>();
    
    return 0;
}
