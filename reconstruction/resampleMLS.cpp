#ifdef _OPENMP
#include <omp.h>
#endif
#include <math/MLS.hpp>
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
#include <boost/math/constants/constants.hpp>
#include <math.h>

#include <format/filename.hpp>
#include <vtk/vtk_utils.hpp>
#include <misc/option_parse.hpp>
#include <image/nrrd_wrapper.hpp>

#include <vtkCellIterator.h>
#include <interpolator.hpp>

#include <dataset.hpp>
#include <mesh.hpp>

#include <Eigen/Eigen>

typedef nvis::bounding_box<nvis::fixed_vector<float, 3> > fbox3;

std::string base(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(0, found);
}

std::string extension(const std::string& name) {
    size_t found = name.find_last_of('.');
    return name.substr(found+1);
}

void load_VTK(celltree::mesh** amesh, celltree::variable** scalars, celltree::variable** vectors, const std::string& filename) {

    VTK_SMART(vtkDataSet) vtk = vtk_utils::readVTK(filename);

    size_t ncells = vtk->GetNumberOfCells();
    celltree::mesh::cell* cells = new celltree::mesh::cell[ncells];

    size_t npts = vtk->GetNumberOfPoints();
    std::vector<float> coords;
    std::vector<unsigned int> ids;

    float* points = new float[npts*3];
    double c[3];
    for (size_t pid=0; pid<npts; ++pid) {
        vtk->GetPoint(pid, c);
        std::copy(c, c+3, &points[3*pid]);
    }

    vtkCellIterator* iter = vtk->NewCellIterator();
    celltree::mesh::cell* cit = cells;
    size_t counter = 0;
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
    {
        vtkIdList *pointIds = iter->GetPointIds();
        int nptsincell = iter->GetNumberOfPoints();
        for (int i=0; i<nptsincell; ++i) {
            ids.push_back(pointIds->GetId(i));
        }
        celltree::cell_kind type;
        switch (iter->GetCellType()) {
            VTK_TETRA: type = celltree::TETRAHEDRON; break;
            VTK_HEXAHEDRON: type = celltree::HEXAHEDRON; break;
            VTK_PYRAMID: type = celltree::PYRAMID; break;
            VTK_WEDGE: type = celltree::PRISM; break;
        }
        (*cit++) = celltree::mesh::cell(type, counter);
        counter += nptsincell;
    }
    unsigned int* indices = new unsigned int[counter];
    std::copy(ids.begin(), ids.end(), indices);

    *amesh = new celltree::mesh(npts, ncells, counter, points, cells, indices);

    *scalars = nullptr;
    if (vtk->GetPointData()->GetScalars() != NULL) {
        float* data = new float[npts];
        for (size_t i=0; i<npts; ++i) {
            data[i] = *vtk->GetPointData()->GetScalars()->GetTuple(i);
        }
        *scalars = new celltree::variable(1, npts, data);
    }

    *vectors = nullptr;
    if (vtk->GetPointData()->GetVectors() != NULL) {
        float* data = new float[3*npts];
        double vec[3];
        for (size_t i=0; i<npts; ++i) {
            vtk->GetPointData()->GetVectors()->GetTuple(i, vec);
            std::copy(vec, vec+3, &data[3*i]);
        }
        *vectors = new celltree::variable(3, npts*3, data);
    }
}

bool load_DLR_or_Nek3D(celltree::mesh** m, celltree::variable** s, celltree::variable** v, const std::string& filename) {
    *m = nullptr;
    *s = nullptr;
    *v = nullptr;
    
    celltree::dataset* ds;
    try {
        ds = celltree::dataset::create(filename);
    }
    catch( std::exception& e) {
        std::cout << "exception caught while importing " << filename << ": " << e.what() << '\n';
        return false;
    }
    *m = ds->read_mesh();
    std::cout << "Mesh read\n";
    *s = ds->read_scalar_variable(0, "pressure");
    std::cout << "Pressure read\n";
    *v = ds->read_vector_variable(0, "velocity");
    std::cout << "Velocity read\n";
    return true;
}

typedef Eigen::Matrix<unsigned int, 3, 1> uint3;
typedef Eigen::Matrix<double, 3, 1> vec3;
typedef Eigen::Matrix<float, 3, 1> fvec3;
typedef Eigen::Matrix<int, 3, 1> int3;
typedef std::pair< vec3, vec3 > box3;
typedef Eigen::Matrix<double, 6, 1> vec6;

// parameters
size_t npts;
uint3 resolution;
double radius;
box3 bounds;
bool verbose, in_parallel;
std::string in_name, out_name, geom_name;

void initialize(int argc, char* argv[])
{
    namespace xcl = spurt::command_line;

    xcl::option_traits
        required(true, false, "Required Options"),
    optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
    "Resample unstructured dataset onto an arbitrary geometry using smooth MLS kernels");
    resolution << 64, 64, 64;
    bounds.first << 0, 0, 0;
    bounds.second << -1, -1, -1;
    verbose = false;
    in_parallel = true;
    vec6 b;
    b << 0, -1, 0, -1, 0, -1;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", in_name, "Input filename", required);
        parser.add_value("output", out_name, "Output base name", required);
        parser.add_value("radius", radius, "MLS radius", required);
        parser.add_tuple<3>("res", resolution, resolution, "Sampling resolution", optional);
        parser.add_value("parallel", in_parallel, in_parallel, "Compute flow map in parallel", optional);
        parser.add_tuple<6>("bounds", b, "Sampling bounds (\"min max\"x3)", optional);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);
        parser.add_value("geometry", geom_name, "Resampling geometry", optional);

        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
            << e.what() << "\n"
                << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    bounds.first << b[0], b[2], b[4];
    bounds.second << b[1], b[3], b[5];
}

template<typename T, int N>
class point {
public:
    typedef T                           scalar_type;
    typedef nvis::fixed_vector<T, N>    vector_type;
    typedef size_t                      index_type;
    typedef T                           value_type;

    point() : m_idx(), m_pos() {}
    point(const index_type& i, const vector_type& v) : m_idx(i), m_pos(v) {}
    point(const point& p) : m_idx(p.idx()), m_pos(p.pos()) {}

    const vector_type& pos() const {
        return m_pos;
    }

    vector_type& pos() {
        return m_pos;
    }
    
    scalar_type* get_pointer() {
        return &m_pos[0];
    }

    const index_type idx() const {
        return m_idx;
    }

    index_type& idx() {
        return m_idx;
    }

    scalar_type distance_to(const point& p) const {
        return norm(m_pos - p.pos());
    }

    scalar_type operator[](size_t i) const {
        return m_pos[i];
    }

private:
    index_type    m_idx;
    vector_type   m_pos;
};

typedef point<double, 3>                  point_type;
typedef point_type::vector_type           vector_type;
typedef KDTree::KDTree<3, point_type>     tree_type;

clock_t init, tree_build_time=0, search_time=0, solve_time=0, total;

fbox3 compute_bounds(const float* points, unsigned int npoints) {
    fbox3 bb;
    const nvis::fvec3* pts = (const nvis::fvec3*)points;
    for (int i=0 ; i<npoints ; ++i) {
        bb.add(pts[i]);
    }
    return bb;
}

template<int NRHS>
void do_fit(const tree_type& tree, const std::vector<double>& rhs, 
            VTK_SMART(vtkDataSet) target, std::vector<double>& soln)
{
    typedef Eigen::Matrix<double, NRHS, 1> value_type;
    
    const value_type* val_rhs = (const value_type*)&rhs[0];

    spurt::MLS::weighted_least_squares<value_type, vector_type> wls(3, 0, NRHS);
    
    soln.resize(NRHS*target->GetNumberOfPoints());
    std::fill(soln.begin(), soln.end(), 0);
    
    int count = 0;
    
    
    std::vector<box3> neighborhoods(target->GetNumberOfPoints());
    std::cout << "Computing size of each vertex's neighborhood\n";
    
    // special case: vtkImageData and descendents: no neighborhood info available.
    if (vtkImageData::SafeDownCast(target) != nullptr) {
        double dx, dy, dz;
        // vtkImageData* = vtkImageData::SafeDownCast(target);
    }
    
    
    for (size_t i=0; i<target->GetNumberOfPoints(); i++) {
        
    }
    

    std::cout << "Fitting in progress" << std::endl;
    for (size_t i=0; i<target->GetNumberOfPoints(); ++i) {
        point_type p;
        target->GetPoint(i, p.get_pointer());
        
        // in radius box find
        std::vector<point_type> in_cube, in_radius;
        clock_t sinit, sfinal;
        sinit=clock();
        tree.find_within_range(p, radius, std::back_inserter(in_cube));
        // std::cout << in_cube.size() << " vertices in cube neighborhood of vertex #" << i << '\n';
        // eliminate points not actually within radius
        for (int j=0; j < in_cube.size(); j++) {
            if (nvis::norm(p.pos()-in_cube[j].pos()) < radius) {
                in_radius.push_back(in_cube[j]);
            }
        }
        if (!in_cube.empty() && in_radius.size() > 5) {
            std::cout << in_radius.size() << " vertices in sphere neighborhood (" 
                    << 100*in_radius.size()/in_cube.size() << "%)\n";
        }
        search_time += clock()-sinit;
        Eigen::MatrixXd fitted_coef;
        std::vector<value_type> rvalues;
        std::vector<vector_type> rpoints;
        if(in_radius.empty()){
            //std::cerr << "No Points within Radius" << std::endl;
            continue;
        }
        // std::cerr << in_cube.size() << " / " << in_radius.size() << " found\n";

        int c_size = in_radius.size();
        rpoints.resize(c_size);
        rvalues.resize(c_size);

        for(int i = 0; i < c_size; i++){
            rpoints[i] = in_radius[i].pos();
            rvalues[i] = val_rhs[in_radius[i].idx()];
        }
        sinit = clock();
        int prec = wls(fitted_coef, rpoints, rvalues, p.pos(), radius);
        solve_time += clock()-sinit;
        for (int n=0; n<fitted_coef.size(); ++n) {
            soln[NRHS*i+n] = fitted_coef(n);
        }
        count++;
    }
}

int main(int argc, char* argv[]) {
    const int dim = 3;

    initialize(argc, argv);
    
    vec3 span = bounds.second - bounds.first;
    bool valid_bounds = std::all_of(&span[0], &span[3], [&](double s) { return s>0; });

    std::cerr << "radius = " << radius << std::endl;

    std::string ext = extension(in_name);

    celltree::mesh* m=nullptr;
    celltree::variable* s=nullptr;
    celltree::variable* v=nullptr;

    if (ext == "dlr" || ext == "nek3d") {
        load_DLR_or_Nek3D(&m, &s, &v, in_name);
    }
    else if (ext == "vtk" || ext == "vtu") {
        load_VTK(&m, &s, &v, in_name);
    }
    else {
        std::cerr << "Unrecognized input filename extension: " << in_name << '\n';
        exit(1);
    }
    
    bool has_scalars = s != nullptr;
    bool has_vectors = v != nullptr;
    bool has_geometry = !geom_name.empty();
    
    VTK_SMART(vtkDataSet) target_mesh;
    if (has_geometry) {
        target_mesh = vtk_utils::readVTK(geom_name);
    }
    else {
        const fbox3 b = compute_bounds(m->points, m->npoints);
        std::cout << "bounds read\n";

        if (!valid_bounds) {
            bounds.first << b.min()[0], b.min()[1], b.min()[2];
            bounds.second << b.max()[0], b.max()[1], b.max()[2];
        }
        std::cout << "bounds = \n" << bounds.first << '\n' << bounds.second << '\n';
        
        vec3 step = ((bounds.second-bounds.first).array() / (resolution.cast<double>() - vec3::Constant(1)).array()).matrix();
        target_mesh = VTK_SMART(vtkUniformGrid)::New();
        vtkUniformGrid::SafeDownCast(target_mesh)->SetDimensions(resolution[0], resolution[1], resolution[2]);
        vtkUniformGrid::SafeDownCast(target_mesh)->SetSpacing(step[0], step[1], step[2]);
        vtkUniformGrid::SafeDownCast(target_mesh)->SetOrigin(bounds.first[0], bounds.first[1], bounds.first[2]);
    }

    std::cout << "dataset: " << m->npoints << " points\n"
              << "         " << m->ncells << " cells\n";

    size_t npts = m->npoints;
    int nrhs = (has_scalars ? 1 : 0) + (has_vectors ? 3 : 0);
    int offset = 0;
    if (has_scalars && has_vectors) offset = 1;
    
    tree_type tree;
    
    clock_t init = clock();

    std::vector<double> values(npts * nrhs);
    for(int i = 0; i < npts; i++){
        double p[3] = { m->points[3*i], m->points[3*i+1], m->points[3*i+2] };
        if (has_scalars) {
            values[nrhs*i] = s->data[i];
        }
        if (has_vectors) {
            std::copy(&v->data[3*i], &v->data[3*i+3], &values[nrhs*i+offset]);
        }
        vector_type pnt(p[0], p[1], p[2]);
        tree.insert(point_type(i, pnt));
    }
    clock_t tree_build_time = clock()-init;

    int n_samples = target_mesh->GetNumberOfPoints();
    std::vector<double> solution;
    
    switch (nrhs) {
        case 1: do_fit<1>(tree, values, target_mesh, solution); break;
        case 3: do_fit<3>(tree, values, target_mesh, solution); break;
        case 4: do_fit<4>(tree, values, target_mesh, solution); break;
        default: {
            std::cerr << "Invalid number of RHS columns: " << nrhs << '\n';
            std::cerr << "Only scalar and/or vector attributes are supported\n";
            exit(1);
        }
        
    }
    clock_t elapsed = clock()-init;
    std::cout << "Algorithm time : " 
        << (double)elapsed / ((double)CLOCKS_PER_SEC) 
        << " s. for " << n_samples << " samples, on "
        << npts << " points"  <<std::endl;
    std::cout << "Tree building time: " 
        << (double)tree_build_time / (double)CLOCKS_PER_SEC << " s.";
    std::cout << search_time << " seconds spent finding points" << std::endl;
    std::cout << solve_time << " seconds spent doing SVD" << std::endl;


    VTK_SMART(vtkDoubleArray) scalars, vectors;
    if (has_scalars) {
        scalars = VTK_SMART(vtkDoubleArray)::New();
        scalars->SetName("pressure");
        scalars->SetNumberOfComponents(1);
        scalars->SetNumberOfTuples(target_mesh->GetNumberOfPoints());
        for (int i=0; i<npts; ++i) {
            scalars->SetTuple1(i, solution[i*nrhs]);
        }
        target_mesh->GetPointData()->SetScalars(scalars);
    }
    if (has_vectors) {
        vectors = VTK_SMART(vtkDoubleArray)::New();
        vectors->SetName("velocity");
        vectors->SetNumberOfComponents(3);
        vectors->SetNumberOfTuples(target_mesh->GetNumberOfPoints());
        for (int i=0; i<npts; ++i) {
            vectors->SetTuple(i, &solution[i*nrhs+offset]);
        }
        target_mesh->GetPointData()->SetVectors(vectors);
    }
    
    vtk_utils::saveVTK_XML(target_mesh, out_name);
    return 0;
}
