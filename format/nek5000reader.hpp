#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>

#include <math/fixed_vector.hpp>
#include <stdexcept>
#include <sstream>

#include "avtNek5000FileFormat.hpp"

namespace spurt {

struct vec_equal {
    bool operator()(const nvis::fvec3& a, const nvis::fvec3& b) const {
        return nvis::all(a == b);
    }
};

class nek5000reader {
    char *metafile, *outfile;
    int timestep;
    bool do_grid;

public:
    nek5000reader(const std::string& metafile) : __fmt(metafile.c_str()) {

    }

    void read_grid(std::vector<unsigned int>& hexind, std::vector<float>& hexpts,
                   std::vector<unsigned int>& rindex) {
        hexind.clear();
        hexpts.clear();

        __fmt.GetMesh(hexpts, hexind);

        int npts = hexpts.size() / 3;
        std::cout << npts << " points in total\n";

        nvis::fvec3 *pts = (nvis::fvec3*) & hexpts.front();

        std::vector<nvis::fvec3> spts(pts, pts + npts);

        std::sort(spts.begin(), spts.end(), nvis::lexicographical_order());
        spts.resize(std::unique(spts.begin(), spts.end(), vec_equal()) - spts.begin());

        std::cout << spts.size() << " unique points = " << (100*spts.size()) / npts << "%%\n";

        std::cout << hexind.size() / 8 << " hexahedra\n";

        std::cout << "Reindexing..." << std::flush;
        std::vector<unsigned int> old2new(npts); // old to new indices
        std::vector<unsigned int> new2old(spts.size()); // new to old indices

        for (unsigned int i = 0; i < npts; ++i) {
            old2new[i] = std::lower_bound(spts.begin(), spts.end(), pts[i], nvis::lexicographical_order()) - spts.begin();
            new2old[old2new[i]] = i;
        }
        for (std::vector<unsigned int>::iterator ii = hexind.begin(); ii != hexind.end(); ++ii)
            *ii = old2new[*ii];
        std::cout << " done\n";

        hexpts = std::vector<float>();
        hexpts.resize(spts.size()*3);
        std::copy(spts.begin(), spts.end(), (nvis::fvec3*)&hexpts[0]);
        rindex.swap(new2old);
        std::cout << rindex.size() << " entries in new2old. Min=" 
            << *std::min_element(rindex.begin(), rindex.end()) << ". Max="
            << *std::max_element(rindex.begin(), rindex.end()) << '\n';
    }

    void read_vectors(std::vector<nvis::fvec3>& vals, int timestep,
                      const std::vector<unsigned int>& rindex) {
        std::vector<float> vec;
        __fmt.GetVectorVar(timestep, vec);

        std::cout << rindex.size() << " vector values to be read\n";

        size_t npts = rindex.size();
        vals.resize(npts);
        std::cout << "reindexing vector field... " << std::flush;
        for (unsigned int i = 0; i < npts; ++i) {
            for (unsigned int j = 0; j < 3; ++j)
                vals[i][j] = vec[3*rindex[i] + j];
        }
        std::cout << "done\n";
    }

private:
    avtNek5000FileFormat __fmt;
};

} // spurt

//---------------------------------------------------------------------------
