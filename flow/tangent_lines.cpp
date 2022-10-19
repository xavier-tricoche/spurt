#include <iostream>
#include <sstream>
#include <ctime>
#include <chrono>

#include <boost/numeric/odeint.hpp>

#include <image/probe.hpp>
#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>
#include <VTK/vtk_utils.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <data/raster.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <format/filename.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

typedef double value_t;

constexpr value_t PI=3.14159265358979323844;
constexpr value_t TWO_PI=6.28318530717958647688;

typedef nvis::fixed_vector< value_t, 3 > state_t;
typedef nvis::fixed_vector< value_t, 12 > ext_state_t;
typedef Eigen::Matrix< value_t, 3, 3 > matrix_t;
typedef Eigen::Matrix< value_t, 3, 1 > col_t;
typedef Eigen::JacobiSVD< matrix_t > svd_t;
typedef nvis::fixed_vector< value_t, 3 > pos_t;
typedef nvis::bounding_box< pos_t > bbox_t;
typedef std::vector< state_t > line_t;
typedef std::vector< value_t > attr_t;

std::string name_in, name_out;
value_t l_max=-1., t_max=10., eps=1.0e-8, dt=1.0e-2, min_dist=1.0e-2;
std::array<size_t, 3> res({ 256, 256, 256 });
std::array< value_t, 6 > bnds({ -PI, PI, -PI, PI, -PI, PI});
bbox_t region;
int nlines=10, fieldid=1;
bool verbose=false, force_2d=false;
bool forward_intg;

enum {
   RANDOM=0,
   X,
   Y,
   Z, 
};

int seeding=Y;


bbox_t to_bbox(const std::array<value_t, 6>& array) {    
    bbox_t b;
    b.min()=pos_t(array[0], array[2], array[4]);
    b.max()=pos_t(array[1], array[3], array[5]);
    return b;
}

template<typename T, size_t N>
nvis::fixed_vector<T, N> to_vec(const std::array<T, N>& array) {
    nvis::fixed_vector<T, N> v;
    for (size_t i=0; i<N; ++i) v[i]=array[i];
    return v;
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map of ABC flow");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required_group);
        parser.add_value("output", name_out, "Output filename", required_group);
        parser.add_value("tmax", t_max, t_max, "Integration length", optional_group);
        parser.add_value("lmax", l_max, l_max, "Maximum curve length", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_value("mind", min_dist, min_dist, "Minimum sampling distance along curves", optional_group);
        parser.add_value("dt", dt, dt, "Temporal sampling", optional_group);
        parser.add_value("number", nlines, nlines, "Number of samples", optional_group);
        parser.add_value("seeding", seeding, seeding, "Sampling pattern", optional_group);
        parser.add_value("field", fieldid, fieldid, "Considered singular vector field", optional_group);
        parser.add_tuple<6>("bounds", bnds, bnds, "Sampling bounds", optional_group);
        parser.add_value("2d", force_2d, force_2d, "Constrain computation to 2d", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

inline value_t cl(const col_t& vals) {
    return (vals(0)-vals(1))/(vals(0)+vals(1)+vals(2));
}

inline value_t cp(const col_t& vals) {
    return static_cast<value_t>(2)*(vals(1)-vals(2))/(vals(0)+vals(1)+vals(2));
}

inline value_t cs(const col_t& vals) {
    return static_cast<value_t>(3)*vals(2)/(vals(0)+vals(1)+vals(2));
}

value_t anisotropy(int i, const col_t& svals) {
    int j=(i+1)%3;
    int k=(i+2)%3;
    value_t dij=(svals(i)-svals(j))*(svals(i)-svals(j));
    value_t dik=(svals(i)-svals(k))*(svals(i)-svals(k));
    return sqrt(dij+dik)/(svals(0)+svals(1)+svals(2));
    
    // switch (i) {
    //     case 0: return cl(svals);
    //     case 1: return cp(svals);
    //     case 2: return cs(svals);
    //     default: throw std::runtime_error("Invalid singular value index");
    // }
}

inline pos_t random_vec() {
    return pos_t(drand48(),
                 drand48(), 
                 drand48());
}

struct NrrdField {
    
    NrrdField(Nrrd* nin) 
        : wrapper(nin) {
        wrapper.use_world();
    }
            
    bool vector(const state_t& x, state_t& v) const {
        return wrapper.value(x, v);
    }
            
    bool jacobian(const state_t& x, matrix_t& J) const {
        return wrapper.jacobian(x, J);
    }
    
    xavier::gage_interface::vector_wrapper wrapper;
};

struct IntegrationMemory {
    enum {
        MAJOR=0,
        MEDIUM,
        MINOR
    };
    
    IntegrationMemory(int ef=MEDIUM) : m_ef(ef), m_valid(false) {}
    
    int m_ef;
    mutable state_t m_last, m_cur;
    mutable value_t m_last_t, m_cur_t;
    mutable bool m_valid;
};

struct NrrdFieldInterface {
    
    NrrdFieldInterface(const NrrdField& field, const IntegrationMemory& memory) 
        : m_field(field), m_memory(memory) {}
    
    NrrdFieldInterface(const NrrdFieldInterface& other) 
        : m_field(other.m_field), m_memory(other.m_memory) {}
    
    void operator()(const state_t& x, state_t& dxdt, value_t t) const {
        matrix_t J;
        bool ok=m_field.jacobian(x, J);
        if (!ok) {
            if (verbose) {
                std::cerr << "WARNING: left domain at " << x << "- unable to interpolate\n";
            }
            std::ostringstream os;
            os << "left domain at " << x;
            throw std::runtime_error(os.str());
        }
        svd_t svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
        value_t aniso=anisotropy(m_memory.m_ef, svd.singularValues());
        col_t singvec=svd.matrixU().col(m_memory.m_ef);
        if (force_2d) singvec[2]=0;
        bool _flipped=false;
        std::ostringstream os;
        
        dxdt=*reinterpret_cast<state_t*>(singvec.data());
        if (m_memory.m_valid) {
            if ((forward_intg && t>m_memory.m_cur_t) || (!forward_intg && t<m_memory.m_cur_t)) {
                m_memory.m_last_t=m_memory.m_cur_t;
                m_memory.m_last=m_memory.m_cur;
                if (verbose)
                    os << "new ref val=" << m_memory.m_cur << " at " 
                        << m_memory.m_cur_t << ". ";
            }
            else {
                if (verbose) 
                    os << "ref val remains " << m_memory.m_last << " at " 
                        << m_memory.m_last_t << ". ";
            }
            if (nvis::inner(dxdt, m_memory.m_last)>0) {
                m_memory.m_cur=dxdt;
                dxdt*=aniso;
                if (verbose)
                    os << " correct orientation. ";
            }
            else {
                _flipped=true;
                m_memory.m_cur=-1.*dxdt;
                dxdt*=-aniso;
                if (verbose)
                    os << " flipped. ";
            }
            m_memory.m_cur_t=t;
        }
        else {
            m_memory.m_valid=true;
            m_memory.m_last=m_memory.m_cur=dxdt;
            m_memory.m_last_t=m_memory.m_cur_t=t;
            dxdt*=aniso;
            if (verbose)
                os << " initialization. ";
        }
        if (verbose) {
            std::cout.setf(std::ios::fixed | std::ios::showpos);
            std::cout << std::setprecision(6)
                      << "\nf(" << x << ", " << t << ")=" 
                      << dxdt << " NB: " << os.str();
        }
    }
    
    const NrrdField& m_field;
    const IntegrationMemory& m_memory;
};

void make_seeds(std::vector<state_t>& seeds) {
    seeds.resize(nlines);
    if (seeding==RANDOM) {
        std::cout << "Random seeding\n";
        srand48(time(0));
        nvis::vec3 sz=region.size();
        sz[2]=0.;
        for (int i=0; i<nlines; ++i) {
            seeds[i]=region.min()+random_vec()*sz;
        }
    }
    else if (seeding==X) {
        std::cout << "X seeding\n";
        pos_t center=region.center();
        value_t span=region.size()[0];
        std::cout << "center=" << center << ", span=" << span << '\n';
        value_t step=span/static_cast<value_t>(nlines+1);
        for (int i=0; i<nlines; ++i) {
            seeds[i]=pos_t(region.min()[0]+(1+i)*step, center[1], center[2]);
        }
    }
    else if (seeding==Y) {
        std::cout << "Y seeding\n";
        pos_t center=region.center();
        value_t span=region.size()[1];
        std::cout << "center=" << center << ", span=" << span << '\n';
        value_t step=span/static_cast<value_t>(nlines+1);
        for (int i=0; i<nlines; ++i) {
            seeds[i]=pos_t(center[0], region.min()[1]+(1+i)*step, center[2]);
        }
    }
    else if (seeding==Z) {
        std::cout << "Z seeding\n";
        pos_t center=region.center();
        value_t span=region.size()[2];
        std::cout << "center=" << center << ", span=" << span << '\n';
        value_t step=span/static_cast<value_t>(nlines+1);
        for (int i=0; i<nlines; ++i) {
            seeds[i]=pos_t(center[0], center[1], region.min()[2]+(1+i)*step);
        }
    }
    else {
        throw std::runtime_error("Unknown seeding pattern");
    }
}

struct length_reached : public std::runtime_error {
    typedef std::runtime_error base_type;

    length_reached(const std::string& what="") : base_type(what) {}
};

struct store_state
{
    std::vector< state_t >& m_states;
    std::vector< value_t >& m_times;
    
    value_t m_length;
    value_t m_maxl;

    store_state(std::vector< state_t >& states, 
                std::vector< value_t >& times,
                const state_t& init_state,
                const value_t& init_time,
                value_t maxl=-1)
    : m_states(states), m_times(times), m_length(0), m_maxl(maxl) {
        m_states.push_back(init_state);
        m_times.push_back(init_time);
    }

    void operator()(const state_t &x, value_t t)
    {
        value_t l=xavier::vector::distance(x, m_states.back());
        if (m_states.empty() || l>min_dist) {
            m_states.push_back(x);
            m_times.push_back(t);
            m_length+=l;
            if (m_maxl>0 && m_length>m_maxl) 
                throw length_reached("requested integration length achieved");
        }
    }
};

std::pair<value_t, value_t> axis_bounds(const NrrdAxisInfo& axis) {
    if (axis.min!=AIR_NAN) 
        return std::make_pair(axis.min, axis.min+axis.size*axis.spacing);
    else if (axis.max!=AIR_NAN) 
        return std::make_pair(axis.max-axis.size*axis.spacing, axis.max);
    else {
        std::cout << "Unable to determine axis bounds\n";
        return std::make_pair(0,0);
    }
}

int main(int argc, const char* argv[])
{
    using namespace xavier;
    using namespace odeint;
    
    initialize(argc, argv);
    
    name_out=xavier::filename::remove_extension(name_out);
    
    region.min()=pos_t(bnds[0], bnds[2], bnds[4]);
    region.max()=pos_t(bnds[1], bnds[3], bnds[5]);
    
    Nrrd* nin=xavier::nrrd_utils::readNrrd(name_in);
    for (int i=0; i<3; ++i) {
        auto min_max=axis_bounds(nin->axis[i+1]);
        if (min_max.first!=min_max.second) {
            region.min()[i]=std::max(region.min()[i], min_max.first);
            region.max()[i]=std::min(region.max()[i], min_max.second);
        }
    }
            
    if (verbose) std::cout << "nb seed points = " << nlines << '\n';
    
    size_t nthreads=1;
    
#if _OPENMP
    if (verbose) std::cout << (nthreads=omp_get_max_threads()) << " threads available\n";
#endif
    
    // storage for tangent lines
    std::vector< line_t > lines(2*nlines);
    std::vector< attr_t > times(2*nlines);
    std::vector< line_t > directions(2*nlines);
    std::vector< attr_t > norms(2*nlines);
    
    std::vector< pos_t > seeds;
    make_seeds(seeds);
    
    int counter = 0;
    int nb_early = 0;
    
    // create a stepper
    runge_kutta_dopri5<state_t> stepper;
    
    if (verbose) std::cout << "bounds=" << region << '\n';
    
    // make sure we have valid min/max info on each spatial axis
    for (int i=0; i<3 ; ++i) {
        if (nin->axis[i+1].min==AIR_NAN && nin->axis[i+1].max==AIR_NAN) {
            if (verbose) 
                std::cout << "Neither min nor max value provided for axis #" << i+1 << '\n';
        }
        else if (nin->axis[i+1].spacing==AIR_NAN) {
            if (verbose)
                std::cout << "No spacing information provided for axis #" << i+1 << '\n';
        }
        else if (nin->axis[i+1].center==nrrdCenterUnknown) {
            if (verbose)
                std::cout << "No center information provided for axis #" << i+1 << '\n';
        }
        
        if (verbose) {
            std::cout << "Axis information available for axis #" << i+1 << '\n'
                << "\tmin=" << nin->axis[i+1].min << '\n'
                << "\tmax=" << nin->axis[i+1].max << '\n'
                << "\tcenter=" << nin->axis[i+1].center << '\n';
            std::pair<value_t, value_t> b=axis_bounds(nin->axis[i+1]);
            std::cout << "\tbounds=" << b.first << ", " << b.second << '\n';
        }
    }
    
    NrrdField field(nin);
    typedef std::shared_ptr<NrrdFieldInterface> field_ptr_t;
    std::vector< field_ptr_t > rhss(nthreads);
    std::vector< IntegrationMemory > rhs_mem(nthreads);
    for (size_t i=0; i<nthreads; ++i) {
        rhss[i]=field_ptr_t(new NrrdFieldInterface(field, rhs_mem[i]));
    }
    
    xavier::ProgressDisplay progress(true);
    
    progress.start(nlines, "Computing tangent lines");
    
    std::clock_t clock_begin=std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (size_t n = 0 ; n < nlines ; ++n) {
        
            #pragma omp atomic
            ++counter;
            
#if _OPENMP
            const int thread=omp_get_thread_num();
#else
            const int thread=0;
#endif
            if (!thread) progress.update(n);
             
            try { 
                forward_intg=true;
                rhs_mem[thread].m_valid=false;
                state_t x = seeds[n];
                if (verbose) {
                    std::cout << "\nstarting forward integration at " << x << '\n';
                }
                integrate_const(
                    make_controlled(eps, eps, stepper), *rhss[thread], x,
                    static_cast<value_t>(0), t_max, dt, 
                    store_state(lines[2*n], times[2*n], x, 0, l_max));
            }
            catch(std::exception& e) {
                if (verbose) {
                    std::cout << "caught exception while integrating forward:"
                        << e.what() << '\n';
                }
                #pragma omp atomic
                ++nb_early;
            }
            try { 
                forward_intg=false;
                state_t x = seeds[n];
                rhs_mem[thread].m_valid=false;
                if (verbose) {
                    std::cout << "\nstarting backward integration at " << x << '\n';
                }
                integrate_const(
                    make_controlled(eps, eps, stepper), *rhss[thread], x,
                    t_max, static_cast<value_t>(0), -dt, 
                    store_state(lines[2*n+1], times[2*n+1], x, t_max, l_max));
            }
            catch(std::exception& e) {
                if (verbose) {
                    std::cout << "caught exception while integrating backward: "
                        << e.what() << '\n';
                }
                #pragma omp atomic
                ++nb_early;
            }
            
            state_t v;
            directions[2*n].resize(lines[2*n].size());
            norms[2*n].resize(lines[2*n].size());
            for (int i=0; i<lines[2*n].size(); ++i) {
                rhss[thread]->operator()(lines[2*n][i], v, times[2*n][i]);
                std::cout << "v(" << lines[2*n][i] << ")=" << v << '\n';
                directions[2*n][i]=v;
                norms[2*n][i]=xavier::vector::norm(v);
            }
            directions[2*n+1].resize(lines[2*n+1].size());
            norms[2*n+1].resize(lines[2*n+1].size());
            for (int i=0; i<lines[2*n+1].size(); ++i) {
                rhss[thread]->operator()(lines[2*n+1][i], v, times[2*n+1][i]);
                std::cout << "v(" << lines[2*n+1][i] << ")=" << v << '\n';
                directions[2*n+1][i]=v;
                norms[2*n+1][i]=xavier::vector::norm(v);
            }
        }
    }
    std::clock_t clock_end=std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
    
    std::cout << "\ntotal wall time as measured by std::chrono was " 
              <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms.";
    std::cout << "\ntotal cpu time was " << 1000.*(clock_end-clock_begin)/CLOCKS_PER_SEC << " ms.";
    std::cout << "\n" << nb_early << " curves left domain (" << 50*(double)nb_early/(double)nlines << "\%)\n";
    
    std::vector< value_t > all_norms;
    std::vector< state_t > all_directions;
    size_t total_n=0;
    for (int i=0; i<2*nlines; ++i) total_n += lines[i].size();
    all_norms.reserve(total_n);
    all_directions.reserve(total_n);
    for (int n=0; n<2*nlines; ++n) {
        const line_t& dirs=directions[n];
        const attr_t& mags=norms[n];
        for (int i=0; i<dirs.size(); ++i) {
            all_norms.push_back(mags[i]);
            std::cout << "adding " << dirs[i] << " to all_directions\n";
            all_directions.push_back(dirs[i]);
        }
    }
    
    std::cout << "there are " << all_directions.size() << " entries in "
        "all_directions\n";
    value_t* test=reinterpret_cast<value_t*>(&all_directions[0]);
    value_t vmax=*std::max_element(test, test+3*all_directions.size());
    std::cout << "vmax=" << vmax << '\n';
    
    vtkSmartPointer<vtkPolyData> polydata=vtk_utils::make_polylines(lines, 1.0e-8);
    vtk_utils::add_scalars(polydata, all_norms);
    vtk_utils::add_vectors(polydata, all_directions);
    vtkDataArray* varray=polydata->GetPointData()->GetVectors();
    int npts=polydata->GetNumberOfPoints();
    int nvecs=polydata->GetPointData()->GetVectors()->GetNumberOfTuples();
    int nvals=polydata->GetPointData()->GetVectors()->GetNumberOfValues();
    std::cout << "# points: " << npts << ", # vecs: " << nvecs << ", # vals: " << nvals << '\n';
    double* v0minmax=varray->GetRange(0);
    double* v1minmax=varray->GetRange(1);
    double* v2minmax=varray->GetRange(2);
    std::cout << "vecs dim 0: " << v0minmax[0] << " -> " << v0minmax[1] << '\n';
    std::cout << "vecs dim 1: " << v1minmax[0] << " -> " << v1minmax[1] << '\n';
    std::cout << "vecs dim 2: " << v2minmax[0] << " -> " << v2minmax[1] << '\n';
    
    vtkSmartPointer<vtkDataSetWriter> writer=vtkSmartPointer<vtkDataSetWriter>::New();
    writer->SetInputData(polydata);
    // writer->SetFileTypeToBinary();
    writer->SetFileName((xavier::filename::remove_extension(name_out)+".vtk").c_str());
    writer->Write();
    
    vtkSmartPointer<vtkPolyData> seedspd=vtk_utils::make_points(seeds);
    writer->SetInputData(seedspd);
    writer->SetFileName((xavier::filename::remove_extension(name_out)+"-seeds.vtk").c_str());
    writer->Write();

    return 0;
}
