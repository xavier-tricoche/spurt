#ifndef __XAVIER_LAVD_HPP__
#define __XAVIER_LAVD_HPP__

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <stdexcept>
#include <locale>
#include <iomanip>

// #include <boost/date_time/posix_time/posix_time.hpp>
// #include <boost/date_time/gregorian/gregorian.hpp>

#include <image/probe.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/log_helper.hpp>
#include <misc/progress.hpp>
#include <vtk/vtk_utils.hpp>
#include <data/raster.hpp>

namespace spurt { namespace lavd {

typedef double value_t;

constexpr double PI =     3.14159265358979323844;
constexpr double TWO_PI = 6.28318530717958647688;
constexpr double SECOND = 1;
constexpr double MINUTE = 60;
constexpr double HOUR =   3600;
constexpr double DAY =    86400;
constexpr double WEEK =   604800;
constexpr double MONTH =  2592000;

// log ostream used 
extern spurt::log::dual_ostream _log_;

typedef nvis::fixed_vector< value_t, 1 > vec1;
typedef nvis::fixed_vector< value_t, 2 > vec2;
typedef nvis::fixed_vector< value_t, 3 > vec3;
typedef nvis::fixed_vector< value_t, 4 > vec4;
typedef nvis::fixed_vector< size_t, 1 > lvec1;
typedef nvis::fixed_vector< size_t, 2 > lvec2;
typedef nvis::fixed_vector< size_t, 3 > lvec3;
typedef nvis::bounding_box< vec2 >     bbox_t;

typedef nvis::fixed_vector< long, 2 > coord_t;
typedef nvis::bounding_box< lvec2 >  coord_range_t;

typedef spurt::image< vec3, 3, value_t, size_t > vector_field_t;
typedef spurt::image< value_t, 3, value_t, size_t > scalar_field_t;
typedef spurt::image< value_t, 1, value_t, size_t > scalar_line_t;

// template< typename T=value_t >
// inline T nrrd_value(const Nrrd* nin, size_t n) {};

template< typename T = value_t >
inline T nrrd_value(const Nrrd* nin, size_t n) {
    return spurt::nrrd_utils::nrrd_data_wrapper<value_t>(nin)[n];
}

template<>
inline vec3 nrrd_value< vec3 >(const Nrrd* nin, size_t n) {
    vec3 v(0);
    v[0] = spurt::nrrd_utils::nrrd_data_wrapper<value_t>(nin)[3*n];
    v[1] = spurt::nrrd_utils::nrrd_data_wrapper<value_t>(nin)[3*n+1];
    return v;
}

inline std::pair<double, double> axis_bounds(const NrrdAxisInfo& axis) {
    if (axis.min!=AIR_NAN) 
        return std::make_pair(axis.min, axis.min+(axis.size-1)*axis.spacing);
    else if (axis.max!=AIR_NAN) 
        return std::make_pair(axis.max-(axis.size-1)*axis.spacing, axis.max);
    else {
        _log_(0) << "Unable to determine axis bounds" << std::endl;
        return std::make_pair(0,0);
    }
}

inline void get_spatial_info(bbox_t& bounds, nvis::vec2& spc, const std::string& filename, int offset=0) {
    // compute bounds of Nrrd file
    Nrrd* nin=spurt::nrrd_utils::readNrrd(filename);
    std::pair<value_t, value_t > minmax;
    minmax = axis_bounds(nin->axis[offset]);
    bounds.min()[0] = minmax.first;
    bounds.max()[0] = minmax.second;
    minmax = axis_bounds(nin->axis[offset+1]);
    bounds.min()[1] = minmax.first;
    bounds.max()[1] = minmax.second;
    spc[0] = nin->axis[offset].spacing;
    spc[1] = nin->axis[offset+1].spacing;
}

void print(const Nrrd* nin) {
    typedef const void* address_t;
    
    _log_(3) << "NRRD: @" << (address_t)(nin) << '\n';
    _log_(3) << "dim=" << nin->dim << '\n';
    _log_(3) << "type=" << nin->type << " (" << spurt::nrrd_utils::nrrd_type_name(nin->type)  << ")" << '\n';
    _log_(3) << "data=" << (address_t)(nin->data) << '\n';
    for (int i=0; i<nin->dim ; ++i) {
        _log_(3) << "axis[" << i << "].size=" << nin->axis[i].size << '\n';
        _log_(3) << "axis[" << i << "].spacing=" << nin->axis[i].spacing << '\n';
        _log_(3) << "axis[" << i << "].min=" << nin->axis[i].min << '\n';
        _log_(3) << "axis[" << i << "].max=" << nin->axis[i].min + nin->axis[i].spacing*(nin->axis[i].size-1) << std::endl;
    }
}

template<size_t N>
struct NrrdScalarField
{
    constexpr static size_t dim = N;
    typedef nvis::fixed_vector< value_t, N > pos_t;
    typedef const void* address_t;
    typedef spurt::gage_interface::scalar_wrapper wrapper_t;

    static vec3 make3d(pos_t p) {
        vec3 q(0);
        for (size_t i=0; i<N; ++i) q[i]=p[i];
        return q;
    }

    static Nrrd* make3d(Nrrd* nin) {
        if (dim==3 || nin->dim==3) {
            _log_(2) << "make3d: did not need to copy\n";
            return nin;
        }
        _log_(1) << "lifting NRRD dataset at " << address_t(nin->data) << " to 3D\n";
        size_t sz[3] = {2, 2, 2};
        value_t spc[3] = {1, 1, 1};
        int center[3] = {nrrdCenterNode, nrrdCenterNode, nrrdCenterNode};
        value_t min[3] = {0, 0, 0};
        value_t* data;
        if (dim==1) {
            sz[0]=nin->axis[0].size;
            spc[0]=nin->axis[0].spacing;
            min[0]=nin->axis[0].min;
            data=(value_t *)calloc(sz[0]*2*2, sizeof(value_t));
            for (size_t i=0; i<sz[0]; ++i) {
                data[i] = data[i+sz[0]] = data[i+2*sz[0]] = data[i+3*sz[0]] =
                    nrrd_value(nin, i); /*((value_t *)nin->data)[i];*/
            }
        }
        if (dim==2) {
            sz[0]=nin->axis[0].size;
            sz[1]=nin->axis[1].size;
            spc[0]=nin->axis[0].spacing;
            spc[1]=nin->axis[1].spacing;
            min[0]=nin->axis[0].min;
            min[1]=nin->axis[1].min;
            data = (value_t *)calloc(sz[0]*sz[1]*2, sizeof(value_t));
            size_t offset=sz[0]*sz[1];
            for (size_t i=0; i<offset; ++i) {
                data[i] = data[i+offset] = nrrd_value(nin, i); 
            }
        }
        Nrrd* nout=nrrdNew();
        if (nrrdWrap_nva(nout, data, spurt::nrrd_utils::nrrd_value_traits_from_type<value_t>::index, 3, sz)) {
            throw std::runtime_error(spurt::nrrd_utils::error_msg("Unable to create Nrrd"));
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, center);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, min);
        return nout;
    }

    NrrdScalarField(Nrrd* nin, const std::string& name="unknown") 
        : m_wrapper(make3d(nin), spurt::gage_interface::BC_INTERP, 
                    false, false, false, false), m_name(name) {
        m_wrapper.use_world();
    #ifdef __VERBOSE_LAVD__
        _log_(2) << "Printing NrrdScalarField " << m_name << ":\n";
        print(nin);
        _log_(2) << std::endl;
    #endif
    }

    ~NrrdScalarField() {}
        
    value_t operator()(const pos_t& x) const {
        value_t s;
    #ifdef __VERBOSE_LAVD__
        _log_(2) << "interpolating " << m_name << " at "
              << x << std::endl;    
    #endif
        if (dim<3) {
            if (!m_wrapper.value(make3d(x), s)) {
                std::ostringstream os;
                os << "Unable to interpolate scalar field at " << x;
                throw std::runtime_error( os.str() );
            }
        }
        else if (!m_wrapper.value(x, s)) {
            std::ostringstream os;
            os << "Unable to interpolate scalar field at " << x;
            throw std::runtime_error( os.str() );
        }
        return s;
    }

    wrapper_t m_wrapper;
    std::string m_name;
};

struct NrrdVectorField
{
    typedef spurt::gage_interface::vector_wrapper wrapper_t;
    typedef wrapper_t::deriv3_t deriv_t;
    
    NrrdVectorField(Nrrd* nin, const std::string name="unknown",
                    bool have_jac=true) 
        : m_wrapper(nin, spurt::gage_interface::BC_INTERP, have_jac),
          m_name(name), m_have_jac(have_jac) {
        m_wrapper.use_world();
    #ifdef __VERBOSE_LAVD__
        _log_(3) << "printing NrrdVectorField " << m_name << ":\n";
        print(nin);
        _log_(3) << std::endl;
    #endif
    }
            
    bool operator()(const vec3& x, vec3& v) const {
        return m_wrapper.value(x, v);
    }
    
    vec3 operator()(const vec3& x) const {
        vec3 v;
        if (!m_wrapper.value(x, v)) {
            std::ostringstream os;
            os << "Probing velocity outside domain of definition at " << x;
            throw std::runtime_error(os.str());
        }
        return v;
    }
    
    value_t vorticity(const vec3& x) const {
        if (!m_have_jac) {
            throw std::runtime_error("Vorticity computation deactivated");
        }
        
        deriv_t J;
        bool ok=m_wrapper.jacobian(x, J);
        if (!ok) {
            std::ostringstream os;
            os << "Probing vorticity outside domain of definition at "
                << x;
            throw std::runtime_error(os.str());
        }
        return J(1,0)-J(0,1);
    }
    
    bool jacobian(const vec3& x, deriv_t& J) const {
        if (!m_have_jac) {
            throw std::runtime_error("Jacobian computation deactivated");
        }
       return m_wrapper.jacobian(x, J);
    }
    
    wrapper_t m_wrapper;
    std::string m_name;
    bool m_have_jac;
};

struct NrrdODERHS
{
    
    typedef NrrdODERHS self_t;
        
    NrrdODERHS( const NrrdVectorField& field, size_t& counter, const bbox_t& region, size_t max_rhs_evals) 
        : m_field(field), m_counter(counter), m_region(region), m_max_rhs_evals(max_rhs_evals) {}
    
    NrrdODERHS( const self_t& other ) 
        : m_field( other.m_field ), m_counter(other.m_counter), m_region(other.m_region),
         m_max_rhs_evals(other.m_max_rhs_evals) {}
    
    void operator()( const vec2& x, vec2& dxdt, value_t t) const {
        if (!m_region.inside(vec2(x[0], x[1]))) {
            std::ostringstream os;
            os << "Interrupting integration because " << x << " is outside selected region: "
                << "min: " << m_region.min() << ", max: " << m_region.max();
            throw std::runtime_error(os.str());
        }
        vec3 v;
        ++m_counter;
        if (m_counter > m_max_rhs_evals) {
            std::ostringstream os;
            os << "\ntoo many RHS queries at " << x;
            throw std::runtime_error(os.str());
        }
        bool ok = m_field( vec3(x[0], x[1], t), v );
        if (!ok) {
            std::ostringstream os;
            os << "\nleft domain at " << x;
            _log_(1) << "WARNING: " << os.str() << std::endl;
            throw std::runtime_error(os.str());
        }
        else if (v[0]==-30000) {
            std::ostringstream os;
            os << "WARNING: reached land mass at " << x;
            throw std::runtime_error(os.str());
        }
        dxdt[0] = v[0];
        dxdt[1] = v[1];
    #ifdef __VERBOSE_LAVD__
        _log_(4) <<  std::setprecision(12) << "v(" <<  x << ", " << t << ")=" << v << std::endl;
    #endif
    }
    
    size_t& m_counter;
    const bbox_t& m_region;
    const NrrdVectorField& m_field;
    const size_t m_max_rhs_evals;
};

size_t compute_support_radius()
{
    _log_(1) << "compute support radius" << std::endl;

    value_t* data = (value_t *)std::calloc(10, sizeof(value_t));
    size_t sz[]={ 10 };
    value_t spc[]={ 1 };
    value_t mins[]={ 0 };
    int center[] = { nrrdCenterNode };
    Nrrd* nin=nrrdNew();
    if (nrrdWrap_nva(nin, data, spurt::nrrd_utils::nrrd_value_traits_from_type<value_t>::index, 
                     1, sz)) {
        throw std::runtime_error(spurt::nrrd_utils::error_msg("unable to compute support radius"));
    }
    nrrdAxisInfoSet_nva(nin, nrrdAxisInfoSpacing, spc);
    nrrdAxisInfoSet_nva(nin, nrrdAxisInfoCenter, center);
    Nrrd* nout = NrrdScalarField<1>::make3d(nin);
    NrrdScalarField<1> sfield(nout, std::string("dummy 1D scalar field in compute_support_radius"));
    size_t r=sfield.m_wrapper.support_radius();
    nrrdNuke(nin);
    return r;
}

inline int to_week(int s) {
    return s/WEEK;
}
inline int to_day(int s) {
    return s/DAY;
}
inline int to_hour(int s) {
    return s/HOUR;
}
inline int to_minute(int s) {
    return s/MINUTE;
}
inline std::string sec2time(int sec) {
    std::ostringstream os;
    int w=to_week(sec);
    int d=to_day(sec-=w*WEEK);
    int h=to_hour(sec-=d*DAY);
    int m=to_minute(sec-=h*HOUR);
    int s=sec-m*MINUTE;
    
    if (w) os << w << " week" << (w>1 ? "s " : " ");
    if (d) os << d << " day" << (d>1 ? "s " : " ");
    if (h) os << h << " hour" << (h>1 ? "s " : " ");
    if (m) os << m << " minute" << (m>1 ? "s " : " ");
    if (s) os << s << " second" << (s>1 ? "s" : "");
    if (os.str().empty()) os << "0";
    
    std::string str = os.str();
    if (str.back()==' ') return str.substr(0, str.size()-1);
    else return str;
}

void validate_value_with_time_unit(value_t& v, const std::string& s, bool allow_times = false) {
    std::string _str;
    size_t n=s.size();
    int i;
    for (i=n-1; i>=0 && std::isalpha(s[i]); --i) {}

    if (i<n-1) {
        _str=s.substr(i+1, std::string::npos);
        if (i>=0) {
            v = std::stod(s.substr(0, i+1));
        }
        else {
            throw std::runtime_error("Invalid time expression: " + s);
        }
    }
    else {
        v=std::stod(s);
        return;
    }

    std::transform(_str.begin(), _str.end(), _str.begin(), ::tolower);
    if (_str=="s" || _str=="sec" || _str=="second" || _str=="seconds") {
    }
    else if (_str=="m" || _str=="min" || _str=="minute" || _str=="minutes") {
        v *= MINUTE;
    } 
    else if (_str=="h" || _str=="hour" || _str=="hours") {
        v *= HOUR;
    }
    else if (_str=="d" || _str=="day" || _str=="days") {
        v *= DAY;
    }
    else if (_str=="w" || _str=="week" || _str=="weeks") {
        v *= WEEK;
    }
    else if (_str=="mo" || _str=="month" || _str=="months") {
        v *= MONTH;
    }
    else if (allow_times && (_str=="x" || _str=="time" || _str=="times")) {
        v *= -1;
    }
    else {
        throw std::runtime_error("unrecognized time unit: " + _str);
    }
}

void validate_date_string(value_t& v, const std::string& s) {
    
    // using namespace boost::posix_time;
    // using namespace boost::gregorian;
    //
    // date earliest(2013,Apr,01);
    // try {
    //     date d(from_simple_string(d));
    //     return
    // }
    // catch(std::exception& e) {
    //
    // }
    
#ifndef __NO_GET_TIME_IN_IOMANIP__    
    std::vector<std::string> date_formats;
    date_formats.push_back("%Y-%b-%d %H:%M:%S");
    date_formats.push_back("%Y-%b-%d %H:%M");
    date_formats.push_back("%Y-%m-%d %H:%M:%S");
    date_formats.push_back("%Y-%m-%d %H:%M");
    date_formats.push_back("%Y-%b-%d");
    date_formats.push_back("%y-%b-%d %H:%M:%S");
    date_formats.push_back("%y-%m-%d %H:%M:%S");
    date_formats.push_back("%y-%b-%d");
    date_formats.push_back("%b %d %Y %H:%M:%S");
    date_formats.push_back("%m/%d %Y %H:%M:%S");
    
    date_formats.push_back("%Y-%b-%d %Hh");
    date_formats.push_back("%Y-%b-%d %Hh");
    date_formats.push_back("%Y-%m-%d %Hh");
    date_formats.push_back("%Y-%m-%d %Hh");
    date_formats.push_back("%Y-%b-%d");
    date_formats.push_back("%y-%b-%d %Hh");
    date_formats.push_back("%y-%m-%d %Hh");
    date_formats.push_back("%y-%b-%d");
    date_formats.push_back("%b %d %Y %Hh");
    date_formats.push_back("%m/%d %Y %Hh");
    
    date_formats.push_back("%b %d %Y");
    date_formats.push_back("%m/%d %Y");
    date_formats.push_back("%b %d %Y %Hh");
    date_formats.push_back("%m/%d %Y %Hh");
    date_formats.push_back("%b %d");
    date_formats.push_back("%m/%d");
    date_formats.push_back("%b %d %Hh");
    date_formats.push_back("%m/%d %Hh");
    
    // try parsing as date in [M]M/[D]D format
    std::tm t;
    t.tm_year = 103;
    t.tm_sec = t.tm_min = t.tm_hour = 0;
    
    for (size_t i=0; i<date_formats.size(); ++i) {
        std::istringstream iss(s);
        iss.imbue(std::locale("en_US"));
        iss >> std::get_time(&t, date_formats[i].c_str());
        if (!iss.fail()) break;
    }
    
    std::tm first;
    first.tm_mday=1;
    first.tm_mon=3;
    first.tm_year=103;
    first.tm_sec = first.tm_min = first.tm_hour = 0;
    std::tm t_copy = t;
    t_copy.tm_isdst=0;
    std::tm first_copy = first;
    first_copy.tm_isdst=0;
    v = std::difftime(std::mktime(&t_copy), std::mktime(&first_copy));
    std::cout << "number of seconds between " << std::asctime(&t) << " and " << std::asctime(&first) << " is "
        << v << std::endl << " that is " << sec2time(v) << '\n';
#else
    throw std::runtime_error("std::get_time not available, dates cannot be parsed");
#endif
}

Nrrd* create_nrrd_volume(std::vector< Nrrd* >& slices, value_t min, value_t spc) {
    Nrrd* nout = nrrdNew();
    _log_(1) << "slices contain " << slices.size() << " nrrds\n";
    size_t sz=1;
    size_t npoints=slices[0]->axis[slices[0]->dim-2].size*slices[0]->axis[slices[0]->dim-1].size;
    bool is_vector_valued=slices[0]->dim==3;
    if (is_vector_valued) {
        sz=3*npoints;
    }
    else {
        sz=npoints;
    }
    value_t* data=(value_t*)calloc(sz*slices.size(), sizeof(value_t));
    for (size_t i=0; i<slices.size(); ++i) {
        if (is_vector_valued) {
            for (size_t n=0; n<npoints; ++n) {
                data[i*sz+3*n] = nrrd_value(slices[i], 2*n);/*((value_t*)(slices[i]->data))[2*n];*/
                data[i*sz+3*n+1] = nrrd_value(slices[i], 2*n+1); /*((value_t*)(slices[i]->data))[2*n+1];*/
                // 3rd component is 0 by default
            }
        }
        else {
            for (size_t n=0; n<npoints; ++n) {
                data[i*sz+n] = nrrd_value(slices[i], n); /*((value_t*)(slices[i]->data))[n];*/
            }
        }
    }
    int dim=slices[0]->dim+1;
    size_t sizes[dim];
    if (is_vector_valued) {
        sizes[0]=3;
        sizes[1]=slices[0]->axis[1].size;
        sizes[2]=slices[0]->axis[2].size;
        sizes[3]=slices.size();
    }
    else {
        sizes[0]=slices[0]->axis[0].size;
        sizes[1]=slices[0]->axis[1].size;
        sizes[2]=slices.size();
    }
    if (nrrdWrap_nva(nout, data, spurt::nrrd_utils::nrrd_value_traits_from_type<value_t>::index, dim, sizes)) {
        throw std::runtime_error("Unable to create NRRD");
    }
    double spacing[dim];
    for (int i=0; i<dim-1; ++i) {
        spacing[i]=slices[0]->axis[i].spacing;
    }
    spacing[dim-1]=spc;
    double mins[dim];
    for (int i=0; i<dim-1; ++i) {
        mins[i]=slices[0]->axis[i].min;
    }
    mins[dim-1]=min;
    int centers[dim];
    for (int i=0; i<dim-1; ++i) {
        centers[i]=nrrdCenterNode;
    }
    centers[dim-1]=nrrdCenterNode;
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spacing);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, mins);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, centers);
    return nout;
}

struct integration_state : public vec4 {
    typedef vec4 base_type;
    
    integration_state(const vec2& p, value_t t, value_t w)
        : base_type(p[0], p[1], t, w) {}
        
        
    vec2& position() { return nvis::subv<0,2,value_t, 4>(*this); }
    const vec2& position() const { return nvis::subv<0,2,value_t, 4>(*this); }
    value_t& time() { return (*this)[2]; }
    value_t time() const { return (*this)[2]; }
    value_t& vorticity() { return (*this)[3]; }
    value_t vorticity() const { return (*this)[3]; }
};

typedef std::vector< integration_state > orbit_type;

struct lavd_state {
    vec2 position;
    value_t time;
    value_t lavd;
    bool stopped;
};

void compute_lavd(std::vector< lavd_state >& lavd_states,
                  const std::vector< orbit_type >& orbits)
{
    size_t N = orbits.size();
    size_t max_nb_steps = 0;
    std::for_each(orbits.begin(), orbits.end(), 
        [&](const orbit_type& orb) {
            max_nb_steps = std::max(max_nb_steps, orb.size());
        });
    
    std::vector< value_t > avg_w(max_nb_steps);
    for (size_t t=0; t<max_nb_steps; ++t) {
        value_t avg = 0;
        size_t n_orbits = 0;
        for (size_t n=0; n<N; ++n) {
            if (orbits[n].size()>t) {
                avg += orbits[n][t].vorticity();
                ++n_orbits;
            }
        }
        avg_w[t] = avg/(value_t)n_orbits;
    }
    
    for (size_t n=0; n<orbits.size(); ++n) {
        if (lavd_states[n].stopped) continue;
        value_t lavd = 0;
        for (size_t t=0; t<orbits[n].size() ; ++t) {
            lavd += fabs(orbits[n][t].vorticity() - avg_w[t]);
        }
        lavd_states[n].time = orbits[n].back().time();
        lavd_states[n].lavd += lavd;
        lavd_states[n].position = orbits[n].back().position();
        lavd_states[n].stopped = orbits.size() == max_nb_steps;
    }
}

template< typename T = float > 
void export_results(value_t current_time, value_t wall_time, value_t cpu_time, 
                    bool final,
                    size_t nb_lost, size_t resx, size_t resy,
                    value_t stepx, value_t stepy,
                    bool export_trajectories,
                    const bbox_t& region,
                    const std::string& file_basename,
                    const std::vector<std::string>& comments,
                    std::vector< lavd_state >& all_states,
                    std::vector< orbit_type >& all_orbits) 
{
    typedef T export_value_t;
    
    size_t nb_samples = resx * resy;
    
    std::string qualifier;
    if (final) qualifier = "Overall";
    else qualifier = "Intermediate";
    
    _log_(1) << "\nEntering export_results\n";
    
    compute_lavd(all_states, all_orbits);

    _log_(0) << "\n" << qualifier << " wall time was " << wall_time << " ms.\n";
    _log_(0) << qualifier << " cpu time was " << cpu_time << " ms.\n";
    _log_(0) << nb_lost << " interrupted trajectories ("
        << static_cast<value_t>(nb_lost)/static_cast<value_t>(nb_samples)*100.
        << "\%)\n";
    
    export_value_t* lavd = (export_value_t *)calloc(4*nb_samples, sizeof(export_value_t));
    _log_(1) << "Filling lavd array in export_results... " << std::flush;
    for (size_t n=0; n<nb_samples; ++n) {
        lavd[4*n  ] = static_cast<export_value_t>(all_orbits[n].back().position()[0]);
        lavd[4*n+1] = static_cast<export_value_t>(all_orbits[n].back().position()[1]);
        lavd[4*n+2] = static_cast<export_value_t>(all_orbits[n].back().time());
        export_value_t val = static_cast<export_value_t>(all_orbits[n].back().vorticity());
        if (std::isnan(val) || std::isinf(val)) {
            val = 0;
        }
        lavd[4*n+3] = val;
    }
    _log_(1) << "done\n";
    
    std::string filename = file_basename + ".nrrd";
    
    _log_(1) << "Setting Nrrd header values... " << std::flush;
    std::vector<size_t> __res(3);
    __res[0] = 4;
    __res[1] = resx;
    __res[2] = resy;
    
    std::vector<double> __mins(3);
    __mins[0] = AIR_NAN;
    __mins[1] = region.min()[0];
    __mins[2] = region.min()[1];
    
    std::vector<double> __spc(3);
    __spc[0] = AIR_NAN;
    __spc[1] = stepx;
    __spc[2] = stepy;
    
    std::vector<int> __ctr(3);
    __ctr[0] = nrrdCenterUnknown;
    __ctr[1] = nrrdCenterNode;
    __ctr[2] = nrrdCenterNode;
        
    _log_(1) << "done\n";
    
    _log_(1) << "Writing NRRD file under " << filename << "... " << std::flush;
    spurt::nrrd_utils::writeNrrdFromContainers(lavd, filename, 
                      __res, __spc, __mins, __ctr, comments);
    _log_(1) << "done\n";
    
    if (export_trajectories) {
        _log_(1) << "Storing trajectories in polydata object... " << std::flush;
        std::vector< std::vector< vec2 > > trajectories(all_orbits.size());
        std::vector<value_t> values;
        for (size_t n=0; n<all_orbits.size(); ++n) {
            trajectories[n].resize(all_orbits[n].size());
            for (size_t t=0; t<all_orbits[n].size(); ++t) {
                trajectories[n][t] = vec2(all_orbits[n][t][0], all_orbits[n][t][1]);
                values.push_back(all_orbits[n][t][3]);
            }
        }
        
        vtkSmartPointer<vtkPolyData> pd(vtk_utils::make_polylines(trajectories, 0.05));
        _log_(1) << "done\n";
        
        _log_(1) << "Adding lavd values to polydata object... " << std::flush;
        vtk_utils::add_scalars(pd, values);
        _log_(1) << "done\n";
        
        vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
        writer->SetInputData(pd);
        filename = file_basename + "_trajectories.vtk";
        _log_(1) << "Writing VTK dataset to file under " << filename << "... " << std::flush;
        writer->SetFileName(filename.c_str());
        writer->SetFileTypeToBinary();
        writer->Write();
        _log_(1) << "done\n";
    }
    
    delete[] lavd;
    
    _log_(1) << "Leaving export_results\n\n";
}

Nrrd* get_border_mask(const std::string& border_mask_name, 
                      const std::string& velocity_name)
{
    if (!border_mask_name.empty()) {
        return spurt::nrrd_utils::readNrrd(border_mask_name);
    }
    
    Nrrd* vel = spurt::nrrd_utils::readNrrd(velocity_name);
    size_t size[2] = { vel->axis[1].size, vel->axis[2].size };
    int *data = (int*)calloc(size[0]*size[1], sizeof(int));
    for (size_t i=0; i<size[0]*size[1]; ++i) {
        double vx = nrrd_value(vel, 2*i);
        double vy = nrrd_value(vel, 2*i+1);
        if (vx!=0 && vy!=0) data[i]=1;
    }
    nrrdNuke(vel);
    
    std::vector<size_t> res(2);
    res[0] = size[0];
    res[1] = size[1];
    
    std::vector<double> mins(2);
    mins[0] = vel->axis[0].min;
    mins[1] = vel->axis[1].min;
    
    std::vector<double> spc(2);
    spc[0] = vel->axis[0].spacing;
    spc[1] = vel->axis[1].spacing;
    
    std::vector<int> ctr(2);
    ctr[0] = nrrdCenterNode;
    ctr[1] = nrrdCenterNode;
    
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, data, nrrdTypeInt, 2, &res[0]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &spc[0]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, &mins[0]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, &ctr[0]);
    
    return nout;
}

void export_mask(const std::string& filename, const bbox_t& region,
                 bool is_scalar=true)
{
    int shift = ( is_scalar ? 0 : 1);

    // export region as mask
    Nrrd* nin(spurt::nrrd_utils::readNrrd(filename));
    value_t minx = nin->axis[shift+0].min;
    value_t miny = nin->axis[shift+1].min;
    value_t dx = nin->axis[shift+0].spacing;
    value_t dy = nin->axis[shift+1].spacing;
  
    size_t mini = (region.min()[0]-minx)/dx;
    size_t maxi = (region.max()[0]-minx)/dx;
    size_t minj = (region.min()[1]-miny)/dy;
    size_t maxj = (region.max()[1]-miny)/dy;
  
    _log_(1) << "index range for selected region: [" << mini << ", " << maxi << "] x [" << minj << ", " << maxj << "]" << std::endl;

    size_t M=nin->axis[shift+0].size;
    size_t N=nin->axis[shift+1].size;
    double* mask_array = (double *)calloc(M*N, sizeof(double));
    for (size_t j=minj; j<=maxj; ++j) {
        for (size_t i=mini; i<maxi; ++i) {
            mask_array[j*M+i]=1;
        }
    }
    size_t sz[2] = {M, N};
    spurt::nrrd_utils::writeNrrd(mask_array, "region_mask.nrrd", nrrdTypeDouble, 2, sz);
    
    nrrdNuke(nin);
    delete[] mask_array;
}


} // lavd
} // spurt


#endif
