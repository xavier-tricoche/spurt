#ifndef __XAVIER_READ_LAMMPS_HPP__
#define __XAVIER_READ_LAMMPS_HPP__

#include <iostream>
#include <fstream>
#include <string>
#include <math/fixed_vector.hpp>
#include <teem/hest.h>
#include <teem/nrrd.h>

namespace xavier {

class LAMMPSreader {
public:
    typedef std::vector<int>    int_vec_t;
    typedef std::vector<float>  float_vec_t;
    typedef std::vector<double> double_vec_t;
    
    LAMMPSreader(const std::string& filename)
        : _filename(filename) {
            
    }
    
private:
    std::string _filename;
    size_t      _timestep;
    size_t      _natoms;
    
    std::vector<std::string> _labels;
    std::vector<std::pair<std::string, int_vec_t> >   int_attributes;
    std::vector<std::pair<std::string, float_vec_t> > float_attributes;
    
};
    
}

char* name_in;
char* base_out;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",    airTypeString,  1,  1,  &name_in,   NULL,   "input file name");
    hestOptAdd(&hopt, "o",      "output",   airTypeString,  1,  1,  &base_out,  NULL,   "output base name");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert LAMMPS dump file to individual NRRD files",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T, int N>
void read(std::fstream& in, T* r)
{
    for (int i=0 ; i<N ; ++i) {
        if (in.eof()) {
            throw std::runtime_error("unexpected end of file");
        }
        in >> r[i];
    }
}

const std::vector< std::string > keywords[] = {
    "TIMESTEP", "NUMBER OF ATOMS", "BOX BOUNDS", "ATOMS"
}

int which_first(std::fstream& in, size_t& pos)
{
    std::string word, first;
    read<std::string, 1>(in, &word);
    for (int i=0 ; i<4 ; ++i) {
        const std::string& kw = keywords[i];
        pos = kw.find(" ");
        first = kw.substr(0, pos);
        if (word == first) {
            return i;
        }
    }
    
    throw std::runtime_error("invalid keyword expression");
    return -1;
}

int find_match(std::fstream& in)
{
    size_t cur, next;
    int id = which_first(in, cur);
    const std::string& kw = keywords[id];
    while (true) {
        next = kw.find(" ", cur+1);
        if (next == std::string::npos) {
            break;
        }
        std::string word = kw.substr(cur, next-cur);
        std::string input;
        read<std::string, 1>(in, &input);
        if (word != input) {
            throw std::runtime_error("invalid keyword expression");
        }
        cur = next;
    }
}

int read_timestep(std::fstream& in, size_t& t)
{
    if (in.eof()) {
        return -1;
    }
    in >> t;
    return 1;
}

int read_nbatoms(std::fstream& in, size_t& n)
{
    if (in.eof()) {
        return -1;
    }
    in >> n;
    return 1;
}

int read_bounds(std::fstream& in, nvis::bounding_box<nvis::fvec3>& box)
{
    float b[6];
    for (int i=0 ; i<6 ; ++i) {
        if (in.eof()) {
            return -1;
        }
        in >> b[i];
    }
    box.min() = nvis::fvec3(b[0], b[2], b[4]);
    box.max() = nvis::fvec3(b[1], b[3], b[5]);
    return 1;
}

int read_atoms(std::fstream& in, float* array, size_t size)
{
    for (size_t i=0 ; i<size ; ++i) {
        if (in.eof()) {
            return -1;
        }
        in >> array[i];
    }
    return 1;
}

int main(int argc, char* argv[])
{

    initialize(argc, argv);
    
    std::fstream in(name_in, std::ios::in);
    if (in.fail()) {
        std::cerr << argv[0] << ": unable to open " << name_in << '\n';
        return -1;
    }
    
    try {
        while (!in.eof()) {
            int ts=-1, natoms=-1;
            nvis::bounding_box<nvis::fvec3> box;
            float* data;
            std::vector<std::string> labels;
            
            while (true) {
                int id = find_match(in);
                switch(id) {
                    case 0:
                        read_timestep(in, ts);
                        break;
                    case 1:
                        read_nbatoms(in, natoms);
                        break;
                    case 2:
                        read_bounds(in, box);
                        break;
                    case 3: {
                        if (ts < 0 || natoms < 0 || labels.empty()) {
                            throw std::runtime_error("missing information");
                        }
                        data = (float*)calloc(labels.size()*natoms, sizeof(float));
                        read_atoms(in, data, labels.size()*natoms);
                        std::ostringstream os;
                    }
                }
            }
        }
    } catch(std::runtime_error& e) {
        std::cerr << "caught: " << e.what() << std::endl;
        return -1;
    }
    
    return 0;
}

#endif
