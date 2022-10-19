#include <iostream>
#include <vector>
#include <string>

#include <math/fixed_vector.hpp>
#include <misc/meta_utils.hpp>
#include <matio.h>

namespace spurt { namespace matlab {
    
template<typename T>
class mat_struct {
public:
    typedef T                   value_type;
    typedef std::vector<T>      array_type;
    typedef std::vector<size_t> coord_type;
    
    mat_struct(const std::string& name=std::string("")) 
        : _name(name), _dims(), _data() {}
        
    mat_struct(const mat_struct& other)
        : _name(other._name), _dims(other._dims), _data(other._data) {}
    
    void resize(const coord_type& dims) {
        _dims=dims;
        data.resize(size(dims));
    }
    
    void set(const std::vector<value_type>& other, const std::vector<size_t> dims) {
        if (other.size() != size(dims)) {
            throw std::runtime_exception("invalid array size");
        }
        _dims=dims;
        _data=other;
    }
    
    array_type& data() { return _data; }
    const array_type data() const { return _data; } 
    
    // meta data lookup
    size_t rank() const { return _dims.size(); }
    const std::string& name() const { return _name; }
    const std::vector<size_t>& dimensions() const { return _dims; }
    bool is_scalar() const { return _data.size()==1; }
    static bool is_complex() const { return std::is_complex<T>::value; }
    static bool is_logical() const { return std::is_same<T, bool>::value; }
    
    // data access
    const value_type& operator[](size_t i, size_t j=0, size_t k=0, size_t l=0) const {
        return _data[idx(i,j,k,l)];
    }
    const value_type& operator[](const coord_type& coord) const {
        return _data[idx(coord)];
    }
    value_type& operator[](size_t i, size_t j=0, size_t k=0, size_t l=0) {
        return _data[idx];
    }
    value_type& operator[](const coord_type& coord) {
        return _data[idx(coord)];
    }
    
private:
    size_t idx(size_t i, size_t j=0, size_t k=0, size_t l=0) const {
        return i + dims[0]*(j + dims[1]*(k + dims[2]*l));
    }

    size_t idx(const coord_type& coord) const {
        size_t r=coord[dimension-1];
        for (int i=dimension-2 ; i>=0 ; --i) {
            r *= dims[i];
            r += coord[i];
        }
        return r;
    }

    static size_t size(const coord_type& d) {
        size_t r=1;
        for (size_t s : dims) r *= s;
        return r;
    }
    
    coord_type              _dims;
    std::vector<value_type> _data;
    std::string             _name;
};

template<typename T> 
std::ostream& operator<<(std::ostream& os, const mat_struct<T>& mat_data) {
    os << "Variable:     " << mat_data.name() << '\n'
       << "\trank:       " << mat_data.rank() << '\n'
       << "\tdimensions: ";
    for (int i=0 ; i<mat_data.rank()-1 ; ++i) {
        os << mat_data.dimensions[i] << "x";
    }
    os << mat_data.dimensions.back() << '\n';
    os << "\tclass:      " << (mat_data.is_scalar() ? "scalar" : "array") << '\n'
       << "\ttype:       "; 
    if (std::is_arithmetic<T>::value) 
        os << type2string<T>::type_name << '\n';
    else 
        os << "non-arithmetic\n";
    os << "\tcomplex?    " << (mat_data.is_complex() ? "yes" : "no") << '\n';
}

mat_t* open_mat(const std::string& name) {
    mat_t *matfp=Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (NULL==matfp) {
        throw std::runtime_error("open_mat: Unable to open " + name);
    }
    return matfp;
}

template<typename T>
void load_mat_info(const std::string& filename, 
                   std::vector<mat_struct<T> >& output) {
    mat_t *matfp=Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (NULL==matfp) {
        throw std::runtime_error("load_mat_info: Unable to open " + name);
    }
    
    matvar_t * matvar;
    while ((matvar=Mat_VarReadNext(matfp)) != NULL) {
        output.push_back(mat_struct(matvar->name));
        mat_struct& _mat=output.back();
    }
    
}

template<typename T>
void load_mat(const std::string& filename, std::vector<mat_struct<T> >& output) {
    mat_t *matfp=Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (NULL==matfp) {
        throw std::runtime_error("load_mat: Unable to open " + filename);
    }
    
    output.clear();
    matvar_t * matvar;
    while ((matvar=Mat_VarReadNext(matfp)) != NULL) {
        output.push_back(mat_struct(matvar->name));
        mat_struct& _mat=output.back();
        
        switch (matvar->data_type) {
            case MAT_T_INT8: 
        }

        if(matvar->data_type==MAT_T_DOUBLE){
            int numel=1;
            for(int i=0; i<matvar->rank; i++){
                numel *= (Vars.back().dim[matvar->rank-1-i]=matvar->dims[i]);
            }
            if(numel==1) Vars.back().type=TYPE_SCALAR;
            else Vars.back().type=TYPE_ARRAY;

            Vars.back().data.resize(numel);
            memcpy(&(Vars.back().data[0]), matvar->data, matvar->nbytes);
        }
        else
            cerr<<"Warning: Type not handled"<<endl;

        
        Mat_VarFree(matvar);
        matvar=NULL;
    }
    
    Mat_Close(matfp);
    
    return true;
}
    
} // matlab
} // xavier
