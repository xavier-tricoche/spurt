#include <crease/measure_wrapper.hpp>
#include <teem/hest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>

struct crease_point_type {
    crease_point_type() {}
    crease_point_type(const nvis::vec3& x) : p(x), val(0), str(0) {}
    
    nvis::vec3 p;
    double val;
    double str;
};

typedef std::list< crease_point_type >          crease_line_type;
typedef crease_line_type::const_iterator        iterator;

void find_best_seed(nvis::vec3& x, double& mode, const crease_line_type& cl,
                    const spurt::MeasureWrapper& measure)
{
    mode = 0;
    for (iterator it = cl.begin() ; it != cl.end() ; ++it) {
        if (fabs(it->val) > mode) {
            mode = fabs(it->val);
            x = it->p;
        }
    }
}

bool correct_point(nvis::vec3& x, const spurt::MeasureWrapper& measure)
{
    double mode = measure.value(x);
    double sign = (mode > 0 ? 1 : -1);
    if (1 - fabs(mode) < 1.0e-6) {
        return true;
    }
    for (int i = 0 ; i < 20 ; ++i) {
        // std::cerr << "mode(" << x << ") = " << mode << std::endl;
        nvis::vec3 dx = measure.gradient(x);
        double m = measure.value(x + dx);
        if (fabs(m) > fabs(mode)) {
            x += dx;
            mode = m;
        } else {
            double h = sign * 0.5;
            while (true) {
                nvis::vec3 y = x + h * dx;
                if (fabs(measure.value(y)) > fabs(mode)) {
                    break;
                }
                h /= 2.;
                if (h < 1.0e-3) {
                    return false;
                }
            }
            x += h * dx;
            mode = measure.value(x);
        }
        
        if (1 - fabs(mode) < 1.0e-6) {
            return true;
        }
    }
    
    return false;
}

double distance(const crease_point_type& cp, const crease_line_type& cl)
{
    double mind = std::numeric_limits<double>::max();
    for (iterator it = cl.begin() ; it != cl.end() ; ++it) {
        double d = nvis::norm(cp.p - it->p);
        if (d < mind) {
            mind = d;
        }
    }
    return mind;
}

double hausdorff(const crease_line_type& cl1, const crease_line_type& cl2)
{
    double max = -std::numeric_limits<double>::max();
    for (iterator it = cl1.begin() ; it != cl1.end() ; ++it) {
        double d = distance(*it, cl2);
        if (d > max) {
            max = d;
        }
    }
    return max;
}

double distance(const crease_line_type& cl, const std::vector<crease_line_type>& cls)
{
    double mind = std::numeric_limits<double>::max();
    for (int i = 0 ; i < cls.size() ; ++i) {
        double d = hausdorff(cl, cls[i]);
        if (d < mind) {
            mind = d;
        }
    }
    return mind;
}


void draw(std::vector<nvis::vec3>& line, const nvis::vec3& seed,
          const spurt::MeasureWrapper& measure, int idx, double h)
{
    const double min_cos = 0.95;
    // std::cerr << "integrating along " << (idx == 0 ? "major" : "minor") << " eigenvector from "
    //           << seed << " where mode is " << measure.value(seed) << std::endl;
    
    std::vector<nvis::vec3> fwd;
    fwd.push_back(seed);
    fwd.push_back(seed + h*measure.eigenvector(seed, idx));
    while (true) {
        nvis::vec3 x = fwd.back();
        nvis::vec3 t = fwd.back() - fwd[fwd.size()-2];
        // std::cout << "after " << fwd.size() << " fwd steps we are at " << x << " where mode is "
        // << measure.value(x) << std::endl;
        t /= nvis::norm(t);
        nvis::vec3 dx = measure.eigenvector(x, idx);
        if (nvis::inner(t, dx) < 0) {
            dx *= -1;
        }
        if (nvis::inner(t, dx) < min_cos) {
            std::cerr << "fwd tracking failed at " << x << " because dot product = " << nvis::inner(t, dx)
                      << std::endl;
            break;
        }
        x += h * dx;
        if (!correct_point(x, measure)) {
            std::cerr << "fwd tracking failed at " << x << " because |1-mode| = " << 1 - fabs(measure.value(x))
                      << std::endl;
            break;
        }
        fwd.push_back(x);
    }
    
    std::vector<nvis::vec3> bwd;
    bwd.push_back(seed);
    bwd.push_back(seed - h*measure.eigenvector(seed, idx));
    while (true) {
        nvis::vec3 x = bwd.back();
        nvis::vec3 t = bwd.back() - bwd[bwd.size()-2];
        // std::cout << "after " << bwd.size() << " bwd steps we are at " << x << " where mode is "
        // << measure.value(x) << std::endl;
        t /= nvis::norm(t);
        nvis::vec3 dx = measure.eigenvector(x, idx);
        if (nvis::inner(t, dx) < 0) {
            dx *= -1;
        }
        if (nvis::inner(t, dx) < min_cos) {
            std::cerr << "bwd tracking failed at " << x << " because dot product = " << nvis::inner(t, dx)
                      << std::endl;
            break;
        }
        x += h * dx;
        if (!correct_point(x, measure)) {
            std::cerr << "bwd tracking failed at " << x << " because mode = " << measure.value(x)
                      << std::endl;
            break;
        }
        bwd.push_back(x);
    }
    
    line.clear();
    for (std::vector<nvis::vec3>::reverse_iterator it = bwd.rbegin() ; it != bwd.rend() ; ++it) {
        line.push_back(*it);
    }
    for (int i = 1 ; i < fwd.size() ; ++i) {
        line.push_back(fwd[i]);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "USAGE: " << argv[0] << " <input.vcl> <output.vcl> <tensor.nrrd>\n";
        return -1;
    }
    
    double h = 0.001;
    if (argc == 5) {
        h = atof(argv[4]);
    }
    
    Nrrd* nrrd = spurt::nrrd_utils::readNrrd(argv[3]);
    spurt::MeasureWrapper measure(nrrd, 2);
    
    std::vector<crease_line_type> all_lines, new_lines;
    
    std::fstream input(argv[1], std::ios::in);
    while (!input.eof()) {
        crease_line_type cl;
        while (!input.eof()) {
            char c;
            float x, y, z;
            input >> c;
            if (c == 'n') {
                break;
            }
            input >> x >> y >> z;
            crease_point_type cp;
            cp.p = nvis::vec3(x, y, z);
            cp.val = measure.value(cp.p);
            cp.str = measure.eigenvalue(cp.p, 1);
            cl.push_back(cp);
        }
        if (cl.size() < 2) {
            continue;
        }
        crease_line_type tmp;
        iterator it = cl.begin();
        nvis::vec3 x, y;
        x = it->p;
        tmp.push_back(*it);
        ++it;
        y = it->p;
        tmp.push_back(*it);
        ++it;
        nvis::vec3 t = y - x;
        t /= nvis::norm(t);
        ++it;
        while (it != cl.end()) {
            x = it->p;
            
            tmp.push_back(*it);
            ++it;
            nvis::vec3 ref = tmp.back().p;
            
            
            all_lines.push_back(cl);
        }
    }
    input.close();
    if (!all_lines.size()) {
        exit(-1);
    }
    bool ridge = all_lines[0].front().val > 0;
    
    for (int i = 0 ; i < all_lines.size() ; ++i) {
        // find best point on current line
        nvis::vec3 seed;
        double mode;
        find_best_seed(seed, mode, all_lines[i], measure);
        if (!correct_point(seed, measure)) {
            continue;
        }
        
        // seed now lies on a degenerate line
        // check if that line has been extracted already
        crease_point_type cp(seed);
        double dist = distance(all_lines[i], new_lines);
        if (dist < 0.01) {
            continue;
        }
        
        // trace line from valid seed point
        std::vector<nvis::vec3> a_line;
        draw(a_line, seed, measure, (ridge ? 0 : 2), h);
        
        // add it to data structure
        crease_line_type new_line;
        for (int i = 0 ; i < a_line.size() ; ++i) {
            new_line.push_back(crease_point_type(a_line[i]));
        }
        new_lines.push_back(new_line);
    }
    
// export results to file
    std::fstream output(argv[2], std::ios::out);
    for (int i = 0 ; i < new_lines.size() ; ++i) {
        for (iterator it = new_lines[i].begin() ; it != new_lines[i].end() ; ++it) {
            output << "p " << it->p[0] << " " << it->p[1] << " " << it->p[2] << '\n';
        }
        output << "n\n";
    }
    output.close();
    
    return 0;
}































