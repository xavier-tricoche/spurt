#ifndef __NRRD_WRAPPER_HPP__
#define __NRRD_WRAPPER_HPP__

#include <teem/nrrd.h>
#include <iostream>
#include <iterator>
#include <memory>
#include <vector>
#include <stdexcept>
#include <array>

#include <math/types.hpp>
#include <math/bounding_box.hpp>
#include <math/vector_manip.hpp>
#include <misc/meta_utils.hpp>
#include <data/raster_data.hpp>


namespace spurt { namespace nrrd_utils {

#define LIST_OF_NRRD_DATA_TYPES \
    X(void,           nrrdTypeUnknown); \
    X(char,           nrrdTypeChar   ); \
    X(unsigned char,  nrrdTypeUChar  ); \
    X(short,          nrrdTypeShort  ); \
    X(unsigned short, nrrdTypeUShort ); \
    X(int,            nrrdTypeInt    ); \
    X(unsigned int,   nrrdTypeUInt   ); \
    X(long,           nrrdTypeLLong  ); \
    X(unsigned long,  nrrdTypeULLong ); \
    X(float,          nrrdTypeFloat  ); \
    X(double,         nrrdTypeDouble );

template<int TypeID>
struct nrrd_value_traits_from_index {};

template<typename T>
struct nrrd_value_traits_from_type {};

#define X(type, nrrdTypeName) \
    template<> \
    struct nrrd_value_traits_from_index<nrrdTypeName> { \
        typedef type data_type; \
        static const int index = nrrdTypeName; \
        static std::string name() { return #nrrdTypeName; } \
    }; \
    template<> \
    struct nrrd_value_traits_from_type<type> { \
        typedef type data_type; \
        static const int index = nrrdTypeName; \
        static std::string name() { return #nrrdTypeName; } \
    };
    LIST_OF_NRRD_DATA_TYPES
#undef X

struct nrrd_value_size {
    size_t operator()(int name) const {
        switch (name) {
        case nrrdTypeChar:
        case nrrdTypeUChar: return 1;
        case nrrdTypeShort:
        case nrrdTypeUShort: return 2;
        case nrrdTypeInt:
        case nrrdTypeUInt: return 4;
        case nrrdTypeLLong:
        case nrrdTypeULLong: return sizeof(long);
        case nrrdTypeFloat: return 4;
        case nrrdTypeDouble: return 8;
        default: throw std::runtime_error("invalid / unknown type");
        }
    }
};

struct nrrd_deleter {
    nrrd_deleter() : m_own(false) {}
    nrrd_deleter(bool own) : m_own(false) {}

    void operator()(Nrrd* nrrd) {
        if (m_own) nrrdNuke(nrrd);
        else nrrdNix(nrrd);
    }

    bool m_own;
};

inline std::string nrrd_type_name(int nrrd_type_index) {
    switch (nrrd_type_index) {
#define X(type, nrrdTypeName) \
    case nrrdTypeName: return #nrrdTypeName;
    LIST_OF_NRRD_DATA_TYPES
#undef X
    default: return "Unknown nrrd type index";
    }
}

inline bool is_volume_data(const Nrrd* nin, bool is_scalar=true) {
    if (is_scalar && nin->dim!=3) return false;
    else if (!is_scalar && nin->dim!=4) return false;
    int offset=(is_scalar ? 0 : 1);
    for (int i=offset; i<offset+3; ++i) {
        if (nin->axis[i].size<2) return false;
    }
    return is_scalar || (!is_scalar && (nin->axis[0].size==3 || nin->axis[0].size==9));
}

template<typename T, size_t N>
class nrrd_params {
public:
    typedef T                          value_type;
    typedef std::array<int, N>         int_array;
    typedef std::array<const char*, N> cstr_array;
    typedef std::array<double, N>      double_array;
    typedef std::array<size_t, N>      long_array;
    typedef std::array<std::string, N> string_array;
    typedef std::vector<std::string>   string_vector;
    static const size_t dim = N;
    static const int data_type = nrrd_value_traits_from_type<T>::index;

    nrrd_params() : m_lbl() {
        m_min.fill(0);
        m_spc.fill(1);
        m_sz.fill(0);
        m_ctr.fill(nrrdCenterUnknown);
    }

    const double_array&  mins() const     { return m_min;  }
    const double_array&  spacings() const { return m_spc;  }
    const long_array&    sizes() const    { return m_sz;   }
    const int_array&     centers() const  { return m_ctr;  }
    const string_array&  labels() const   { return m_lbl;  }
    const string_vector& comments() const { return m_cmts; }

    double_array&  mins()     { return m_min; }
    double_array&  spacings() { return m_spc; }
    long_array&    sizes()    { return m_sz;  }
    int_array&     centers()  { return m_ctr; }
    string_array&  labels()   { return m_lbl; }
    string_vector& comments() { return m_cmts; }

    // C-style access
    const double* C_mins() const {
        return m_min.data();
    }
    const double* C_spacings() const {
        return m_spc.data();
    }
    const size_t* C_sizes() const {
        return m_sz.data();
    }
    const char** C_labels() const {
        static cstr_array clbl;
        for (size_t i=0 ; i<m_lbl.size() ; ++i) {
            clbl[i] = m_lbl[i].c_str();
        }
        return clbl.data();
    }
    const int* C_centers() const {
        return m_ctr.data();
    }

private:
    double_array  m_min, m_spc;
    long_array    m_sz;
    int_array     m_ctr;
    string_array  m_lbl;
    string_vector m_cmts;
};

template<typename T>
inline bool invalid(T v)
{
    return (std::isnan(v) || std::isinf(v));
}

inline size_t nrrd_raster_size(const Nrrd* A, bool is_scalar=true) {
    size_t s=1;
    for (int i=(is_scalar ? 0 : 1); i<A->dim; ++i) {
        s*=A->axis[i].size;
    }
    return s;
}

inline bool matching_raster_sizes(const Nrrd* A, const Nrrd* B, bool is_scalar=true) {
    if (A->dim!=B->dim) return false;
    for (int d=(is_scalar ? 0 : 1); d<A->dim; ++d) {
        if (A->axis[d].size!=B->axis[d].size) return false;
    }
    return true;
}

class nrrd_traits {

    typedef nrrd_traits self_type;

public:
    nrrd_traits() {}

    nrrd_traits(const Nrrd* nrrd) {
        m_dim = nrrd->dim;
        m_type = nrrd->type;
        m_sizes.resize(m_dim);
        m_mins.resize(m_dim, 0);
        m_maxs.resize(m_dim, 0);
        m_spacings.resize(m_dim, 1);
        m_centers.resize(m_dim, nrrdCenterNode);
        m_labels.resize(m_dim);

        if (nrrd->cmtArr->len) {
            m_comments.resize(nrrd->cmtArr->len);
            for (int i=0; i<m_comments.size(); ++i) {
                m_comments[i] = nrrd->cmt[i];
            }
        }

        for (int d=0; d<m_dim; ++d) {
            NrrdAxisInfo axis = nrrd->axis[d];
            m_sizes[d] = axis.size;
            if (!invalid(axis.min)) m_mins[d] = axis.min;
            else if (!invalid(axis.max)) m_maxs[d] = axis.max;
            if (axis.center != nrrdCenterUnknown) m_centers[d] = axis.center;
            if (!invalid(axis.spacing)) m_spacings[d] = axis.spacing;
            if (axis.label) m_labels[d] = axis.label;
            if (m_maxs[d] == 0) m_maxs[d] = m_mins[d] + m_spacings[d]*(m_sizes[d]-1);
        }
    }

    size_t dim() const { return m_dim; }
    int type() const { return m_type; }

    size_t mem_size() const {
        nrrd_value_size vs;
        size_t value_size = vs(m_type);
        size_t nvalues=1;
        std::for_each(m_sizes.begin(), m_sizes.end(), [&](size_t s){
            nvalues *= s;
        });
        return value_size * nvalues;
    }

    const std::vector<double>&      mins() const     { return m_mins;  }
    const std::vector<double>&      maxs() const     { return m_maxs;  }
    const std::vector<double>&      spacings() const { return m_spacings;  }
    const std::vector<size_t>&      sizes() const    { return m_sizes;   }
    const std::vector<int>&         centers() const  { return m_centers;  }
    const std::vector<std::string>& labels() const   { return m_labels;  }
    const std::vector<std::string>& comments() const { return m_comments; }

    std::vector<double>&      mins()      { return m_mins;  }
    std::vector<double>&      maxs()      { return m_maxs;  }
    std::vector<double>&      spacings()  { return m_spacings;  }
    std::vector<size_t>&      sizes()     { return m_sizes;   }
    std::vector<int>&         centers()   { return m_centers;  }
    std::vector<std::string>& labels()    { return m_labels;  }
    std::vector<std::string>& comments()  { return m_comments; }

private:
    size_t m_dim;
    int m_type;
    std::vector<double> m_mins, m_maxs, m_spacings;
    std::vector<size_t> m_sizes;
    std::vector<int> m_centers;
    std::vector<std::string> m_labels;
    std::vector<std::string> m_comments;
};

class nrrd_wrapper {
public:
    nrrd_wrapper() : m_nrrd(NULL), m_own(false ){}
    nrrd_wrapper(Nrrd* nin, bool own=false)
        : m_nrrd(nin, nrrd_deleter(own)), m_own(own) {}
    nrrd_wrapper(const std::string& filename, bool own=true);
    nrrd_wrapper(const nrrd_wrapper& other)
        : m_nrrd(other.m_nrrd), m_own(other.m_own) {}
    ~nrrd_wrapper() {}

    Nrrd* pointer() { return m_nrrd.get(); }
    const Nrrd* pointer() const { return m_nrrd.get(); }

    bool save(const std::string& filename) const {
        return !nrrdSave(filename.c_str(), m_nrrd.get(), NULL);
    }

    nrrd_traits traits() const {
        return nrrd_traits(m_nrrd.get());
    }

private:
    mutable bool m_own;
    std::shared_ptr<Nrrd> m_nrrd;
    nrrd_traits m_traits;
};

inline std::string error_msg(const std::string& fun_name="", const char* what=NRRD) {
    return fun_name + ": " + biffGetDone(what);
}

template<int N>
inline void
compute_raster_bounds(spurt::bounding_box<small_vector<double, N>>& bounds,
                      const Nrrd* nrrd, bool is_scalar=true)
{
    int skip = (is_scalar ? 0 : 1);
    bounds.min() = 0;
    bounds.max() = 0;
    assert(nrrd->dim - skip == N);
    for (int i = skip ; i < nrrd->dim ; ++i) {
        int j = i-skip;
        const NrrdAxisInfo& axis = nrrd->axis[i];
        double step = axis.spacing;
        if (invalid(step)) step = 1;
        double width = (axis.size - 1) * step;
        if (!invalid(axis.min)) {
            bounds.min()[j] = axis.min;
            bounds.max()[j] = axis.min + width;
        }
        else if (!invalid(axis.max)) {
            bounds.max()[j] = axis.max;
            bounds.min()[j] = axis.max - width;
        }
        else {
            bounds.min()[j] = 0;
            bounds.max()[j] = width;
        }
        if (axis.center == nrrdCenterCell) {
            bounds.min()[j] -= 0.5 * step;
            bounds.max()[j] += 0.5 * step;
        }
    }
}

template<int N> inline
spurt::bounding_box<small_vector<double, N>>
get_bounds(const Nrrd* nrrd, bool is_scalar=true)
{
    bounding_box<small_vector<double, N>> b;

    compute_raster_bounds<N>(b, nrrd, is_scalar);
    return b;
}

template<int N> inline
spurt::bounding_box<small_vector<double, N>> get_bounds(const nrrd_wrapper& wrap) {
    return get_bounds<N>(wrap.pointer());
}

template<int N> inline
Eigen::Vector<double, N> step(const Nrrd* nrrd)
{
    Eigen::Vector<double, N> s;
    for (int i = 0 ; i < N ; ++i)
    {
        s[i] = nrrd->axis[nrrd->dim-N+i].spacing;
        if (invalid(s[i])) s[i] = 1;
    }
    return s;
}

template<int N> inline
spurt::bounding_box<small_vector<double, N> >step(const nrrd_wrapper& wrap) {
    return step<N>(wrap.pointer());
}

template<typename T>
struct nrrd_data_wrapper {
    nrrd_data_wrapper(const Nrrd* nrrd) : m_nrrd(nrrd) {}
    nrrd_data_wrapper(const nrrd_wrapper& wrap) : m_nrrd(wrap.pointer()) {}
    T operator[](size_t i) const {
        switch (m_nrrd->type) {
        case nrrdTypeChar:
            return static_cast<T>(((char*)m_nrrd->data)[i]);
        case nrrdTypeUChar:
            return static_cast<T>(((unsigned char*)m_nrrd->data)[i]);
        case nrrdTypeShort:
            return static_cast<T>(((short*)m_nrrd->data)[i]);
        case nrrdTypeUShort:
            return static_cast<T>(((unsigned short*)m_nrrd->data)[i]);
        case nrrdTypeInt:
            return static_cast<T>(((int*)m_nrrd->data)[i]);
        case nrrdTypeUInt:
            return static_cast<T>(((unsigned int*)m_nrrd->data)[i]);
        case nrrdTypeLLong:
            return static_cast<T>(((long int*)m_nrrd->data)[i]);
        case nrrdTypeULLong:
            return static_cast<T>(((unsigned long int*)m_nrrd->data)[i]);
        case nrrdTypeFloat:
            return static_cast<T>(((float*)m_nrrd->data)[i]);
        case nrrdTypeDouble:
            return static_cast<T>(((double*)m_nrrd->data)[i]);
        default:
            throw std::runtime_error("unrecognized data type\n");
        }
    }
    T operator()(size_t i) const {
        return (*this)[i];
    }

    const Nrrd* m_nrrd;
};

template<typename T1, typename T2>
inline void to_vector(std::vector<T1>& vals, const void* data, size_t size)
{
    typedef data_traits<T1> traits_type;
    typedef typename traits_type::value_type value_type;
    size_t nvals = traits_type::size();

    vals.resize(size);
    value_type* ptr = (value_type *)&(vals[0]);
    for (size_t i = 0 ; i < size*nvals ; ++i) {
        (*ptr++) = ((T2*)data)[i];
    }
}

template<typename T>
inline void to_vector(std::vector<T>& vals, const Nrrd* nin)
{
    size_t size = 1;
    for (size_t i = 0 ; i < nin->dim ; ++i) {
        size *= nin->axis[i].size;
    }
    // std::cerr << "size = " << size << std::endl;

    vals.clear();
    switch (nin->type) {
    case nrrdTypeChar:
        return to_vector<T, char>(vals, nin->data, size);
    case nrrdTypeUChar:
        return to_vector<T, unsigned char>(vals, nin->data, size);
    case nrrdTypeShort:
        return to_vector<T, short>(vals, nin->data, size);
    case nrrdTypeUShort:
        return to_vector<T, unsigned short>(vals, nin->data, size);
    case nrrdTypeInt:
        return to_vector<T, int>(vals, nin->data, size);
    case nrrdTypeUInt:
        return to_vector<T, unsigned int>(vals, nin->data, size);
    case nrrdTypeLLong:
        return to_vector<T, long int>(vals, nin->data, size);
    case nrrdTypeULLong:
        return to_vector<T, unsigned long int>(vals, nin->data, size);
    case nrrdTypeFloat:
        return to_vector<T, float>(vals, nin->data, size);
    case nrrdTypeDouble:
        return to_vector<T, double>(vals, nin->data, size);
    default:
        throw std::runtime_error("unrecognized data type\n");
    }
}

template<typename T>
inline void to_vector(std::vector<T>& vals, const nrrd_wrapper& wrap) {
    to_vector<T>(vals, wrap.pointer());
}

template<typename Size_, typename Scalar_, size_t Dim, typename Value_, 
         typename Coord_ = small_vector<Size_, Dim>, typename Pos_ = small_vector<Scalar_, Dim>>
inline void to_raster(raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_>& out, 
                      const Nrrd* nin, bool is_scalar=true) 
{
    typedef raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_> raster_type;
    typedef typename raster_type::grid_type grid_type;
    typedef typename raster_type::coord_type coord_type;
    typedef typename raster_type::value_type value_type;
    
    auto bounds = get_bounds<Dim>(nin, is_scalar);
    nrrd_traits traits(nin);
    int offset = is_scalar ? 0 : 1;
    assert(nin->dim == Dim+offset);
    assert( !is_scalar || std::is_scalar<value_type>::value );
    coord_type res = 0;
    auto _res = traits.sizes();
    for (int i=0; i<Dim; ++i) res[i] = _res[i];
    grid_type agrid(res, bounds);
    out = raster_type(agrid);
    to_vector(out.data(), nin);
}

template<typename T1, typename T2>
inline T1* to_array(const void* data, size_t size)
{
    T1* array = (T1*)calloc(size, sizeof(T1));
    for (int i = 0 ; i < size ; ++i) {
        array[i] = ((T2*)data)[i];
    }
    return array;
}

template<typename T>
inline T* to_array(const Nrrd* nin)
{
    size_t size = 1;
    for (int i = 0 ; i < nin->dim ; ++i) {
        size *= nin->axis[i].size;
    }

    switch (nin->type) {
    case nrrdTypeChar:
        return to_array<T, char>(nin->data, size);
    case nrrdTypeUChar:
        return to_array<T, unsigned char>(nin->data, size);
    case nrrdTypeShort:
        return to_array<T, short>(nin->data, size);
    case nrrdTypeUShort:
        return to_array<T, unsigned short>(nin->data, size);
    case nrrdTypeInt:
        return to_array<T, int>(nin->data, size);
    case nrrdTypeUInt:
        return to_array<T, unsigned int>(nin->data, size);
    case nrrdTypeLLong:
        return to_array<T, long int>(nin->data, size);
    case nrrdTypeULLong:
        return to_array<T, unsigned long int>(nin->data, size);
    case nrrdTypeFloat:
        return to_array<T, float>(nin->data, size);
    case nrrdTypeDouble:
        return to_array<T, double>(nin->data, size);
    default:
        throw std::runtime_error("unrecognized data type\n");
    }
}

template<typename T>
inline T* to_array(const nrrd_wrapper& wrap) {
    return to_array<T>(wrap.pointer());
}

template< typename T >
inline Nrrd* make3d(Nrrd* nin, bool is_scalar=true) {
    typedef T value_t;

    int offset = is_scalar ? 0 : 1;
    if (nin->dim == 3 + offset) return nin;

    size_t sz[4] = {2, 2, 2, 2};
    double spc[4] = {1, 1, 1, 1};
    int ctr[4] = {nrrdCenterNode, nrrdCenterNode, nrrdCenterNode, nrrdCenterNode};
    double min[4] = {0, 0, 0, 0};
    value_t* data;

    nrrd_data_wrapper<value_t> value_wrapper(nin);

    assert(nin->dim == 2 + offset);
    for (int i=0; i<nin->dim; ++i) {
        sz[i] = nin->axis[i].size;
        spc[i] = nin->axis[i].spacing;
        min[i] = nin->axis[i].min;
        ctr[i] = nin->axis[i].center;
    }
    size_t nscalars = is_scalar ? sz[0] * sz[1] : sz[0] * sz[1] * sz[2];
    data = (value_t *)calloc(nscalars * 2, sizeof(value_t));
    size_t stride = nscalars;
    for (size_t i=0; i<stride; ++i) {
            data[i] = data[i+stride] = value_wrapper[i];
    }

    Nrrd* nout=nrrdNew();
    if (nrrdWrap_nva(nout, data, nrrd_value_traits_from_type<value_t>::index, 3+offset, sz)) {
        throw std::runtime_error(error_msg("Unable to create Nrrd"));
    }
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, ctr);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, min);

    return nout;
}


inline Nrrd* readNrrd(const std::string& filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        throw std::runtime_error(error_msg("readNrrd"));
    }
    return nin;
}

inline nrrd_wrapper::nrrd_wrapper(const std::string& filename, bool own) {
    m_nrrd.reset(readNrrd(filename), nrrd_deleter(own));
    m_own = own;
}

inline void writeNrrd(Nrrd* nout, const std::string& filename,
                      bool compressed=false) {
    NrrdIoState *nio = nullptr;
    if (compressed) {
        nio = nrrdIoStateNew();
        nio->encoding = nrrdEncodingBzip2;
        nio->bzip2BlockSize = -1;
        nio->format = nrrdFormatNRRD;
    }
    if (nrrdSave(filename.c_str(), nout, nio)) {
        throw std::runtime_error(error_msg("writeNrrd: error while saving"));
    }
}

inline void writeNrrd(void* data, const std::string& filename,
                      int data_type, int ndim, size_t* dims,
                      bool compressed=false)
{
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, data, data_type, ndim, dims)) {
        throw std::runtime_error(error_msg("writeNrrd: error while wrapping"));
    }
    writeNrrd(nout, filename, compressed);
    nrrdNix(nout); // only deletes the structure not the data
}

struct comment_inserter {
    static Nrrd* add(Nrrd* nin, const std::string& comment) {
        nrrdCommentAdd(nin, comment.c_str());
        return nin;
    }
    static Nrrd* add(Nrrd* nin, const std::vector<std::string>& comments) {
        std::for_each(comments.begin(), comments.end(),
            [&](const std::string& comment) {
                nrrdCommentAdd(nin, comment.c_str());
            }
        );
        return nin;
    }
};

template< typename DataType_, typename Container1_,
          typename Container2_=std::vector<double>,
          typename Container3_=std::vector<double>,
          typename Container4_=std::vector<int>,
          typename Comments_=std::string >
inline void writeNrrdFromContainers(DataType_* data,
                                    const std::string& filename,
                                    const Container1_& dims,
                                    const Container2_& spc = Container2_(),
                                    const Container3_& min = Container3_(),
                                    const Container4_& center = Container4_(),
                                    const Comments_& comments = Comments_())
{
    Nrrd *nout = nrrdNew();
    int index = nrrd_value_traits_from_type<DataType_>::index;

    std::vector< size_t > __dims(dims.begin(), dims.end());
    std::vector< double > __spc(dims.size());
    std::vector< double > __min(dims.size());
    std::vector< int >    __ctr(dims.size());

    if (nrrdWrap_nva(nout, data, index, dims.size(), &__dims.front())) {
        throw std::runtime_error(error_msg("writeNrrd"));
    }

    if ( spc.size() ) {
        std::copy(spc.begin(), spc.end(), __spc.begin());
    }
    else {
        std::fill(__spc.begin(), __spc.end(), 1.);
    }
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc.front());

    if ( min.size() ) {
        std::copy(min.begin(), min.end(), __min.begin());
    }
    else {
        std::fill(__min.begin(), __min.end(), 0.);
    }
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, &__min.front());

    if ( center.size() ) {
        std::copy(center.begin(), center.end(), __ctr.begin());
    }
    else {
        std::fill(__ctr.begin(), __ctr.end(), nrrdCenterNode);
    }
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, &__ctr.front());

    comment_inserter::add(nout, comments);

    if (nrrdSave(filename.c_str(), nout, NULL)) {
        throw std::runtime_error(error_msg("writeNrrd"));
    }

    nrrdNix(nout);  // only deletes the structure not the data
}

template<typename T, size_t N>
inline void writeNrrdFromParams(void* data, const std::string& filename,
                                const nrrd_params<T, N>& params)
{
    size_t dim = nrrd_params<T, N>::dim;

    Nrrd *nout = nrrdNew();

    if (nrrdWrap_nva(nout, data, params.data_type, dim, params.C_sizes())) {
        throw std::runtime_error(error_msg("writeNrrd"));
    }
    // Per-axis info
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, params.C_spacings());
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, params.C_mins());
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, params.C_centers());
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, params.C_labels());
    // Comments to store ad-hoc meta-data
    for (size_t i=0 ; i<params.comments().size() ; ++i) {
        nrrdCommentAdd(nout, params.comments()[i].c_str());
    }
    if (nrrdSave(filename.c_str(), nout, NULL)) {
         throw std::runtime_error(error_msg("writeNrrd"));
    }

    nrrdNix(nout);  // only deletes the structure not the data
}

} // nrrd_utils

} // spurt

#endif
