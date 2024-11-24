#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <math/types.hpp>
#include <math/bounding_box.hpp>
#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <misc/progress.hpp>

#include <format/filename.hpp>
#include <vtk/vtk_utils.hpp>

// celltree
#include <celltree.hpp>
#include <celltree_builder.hpp>
#include <dataset.hpp>
#include <interpolator.hpp>
#include <mesh.hpp>

#include <vtkCellIterator.h>

typedef spurt::bounding_box<spurt::fvec3> fbox3;

// parameters
spurt::ivec3   res;
std::string   in_name, out_name;
fbox3         bounds;

std::string me;
void usage( const std::string& msg="")
{
    if (!msg.empty()) {
        std::cerr << "ERROR: " << msg << std::endl;
    }
    std::cout
        << "Usage  : " << me << " [parameters] [options]\n"
        << '\n'
        << "Synopsis: resample given dataset using cell tree and native interpolation"
        << '\n'
        << "Parameters:\n"
        << " -i  | --input <string>       Dataset info file\n"
        << " -o  | --output <string>      Output file name\n"
        << " -r  | --resolution <int> x3  Resampling resolution in X, Y, and Z\n"
        << " -g  | --grid <string>        Use provided resampling grid\n"
        << '\n'
        << "Options:\n"
        << " -b  | --bounds <float> x6    Bounding box of region (in world coordinates).\n"
        << "                              Default: all\n"
        << std::endl;

    exit(!msg.empty());
}

fbox3 compute_bounds(const float* points, unsigned int npoints) {
    fbox3 bb;
    const spurt::fvec3* pts = (const spurt::fvec3*)points;
    for (int i=0 ; i<npoints ; ++i) {
        bb.add(pts[i]);
    }
    return bb;
}

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
        // vtkPoints *points = iter->GetPoints();
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

int main(int argc, char* argv[]) {
    typedef celltree::interpolator::coord_type coord_type;
    typedef celltree::interpolator::value_type value_type;

    res = spurt::ivec3(0);
    std::string grid_name;

    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                usage("missing input");
            }
            in_name = argv[++i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                usage("missing output");
            }
            out_name = argv[++i];
        }
        else if (arg == "-h" || arg == "--help") {
            usage();
        }
        else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-3) {
                usage("missing resolution information");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
            res[2] = atoi(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-6) {
                usage("missing bounds information");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            bounds.min()[2] = atof(argv[++i]);
            bounds.max()[2] = atof(argv[++i]);
        }
        else if (arg == "-g" || arg == "--grid") {
            if (i == argc-1) {
                usage("missing grid filename");
            }
            grid_name = argv[++i];
        }
        else {
            usage("invalid argument: " + arg);
        }
    }

    if (in_name.empty() || out_name.empty()) {
        usage("missing input or output file name");
    }
    else if (*std::min_element(res.begin(), res.end()) <= 0 && grid_name.empty()) {
        usage("missing / invalid resolution information");
    }

    // user reader appropriate for input file type
    spurt::timer _timer;

    std::string ext = extension(in_name);

    celltree::mesh* m=nullptr;
    celltree::variable* s=nullptr;
    celltree::variable* v=nullptr;

    if (ext == "dlr" || ext == "nek3d") {
        celltree::dataset* ds;
        try {
            ds = celltree::dataset::create(in_name);
        }
        catch( std::exception& e) {
            std::cout << "exception caught while importing " << in_name << ": " << e.what() << '\n';
        }
        m = ds->read_mesh();
        std::cout << "Mesh read\n";
        s = ds->read_scalar_variable(0, "pressure");
        std::cout << "Pressure read\n";
        v = ds->read_vector_variable(0, "velocity");
        std::cout << "Velocity read\n";
    }
    else if (ext == "vtk" || ext == "vtu") {
        load_VTK(&m, &s, &v, in_name);
    }
    else {
        std::cerr << "Unrecognized input filename extension: " << in_name << '\n';
        exit(1);
    }

    spurt::fvec3 step;

    bool has_geometry = !grid_name.empty();
    VTK_SMART(vtkDataSet) target_mesh;
    if (has_geometry) {
        target_mesh = vtk_utils::readVTK(grid_name);
    }
    else {
        const fbox3 b = compute_bounds(m->points, m->npoints);
        std::cout << "bounds read\n";

        if (spurt::any(bounds.min() >= bounds.max())) bounds = b;
        std::cout << "bounds = " << bounds << '\n';

        step = bounds.size() / (res - 1);
    }

    std::cout << "dataset: " << m->npoints << " points\n"
              << "         " << m->ncells << " cells\n";

    spurt::timer timer;
    celltree::celltree ct;
    {
        celltree::celltree_builder builder;
        builder.m_buckets  = 5;
        builder.m_leafsize = 8;
        builder.build( ct, *m );
    }
    std::cout << "celltree built in " << timer.elapsed() << " seconds\n";
    celltree::interpolator* vintp;
    if (v != nullptr) {
        vintp = new celltree::interpolator( m, v, ct );
    }
    celltree::interpolator* sintp;
    if (s != nullptr) {
        sintp = new celltree::interpolator( m, s, ct );
    }

    size_t number_of_threads = 1;
#ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
#endif

    size_t number_of_samples;
    if (has_geometry) {
        number_of_samples = target_mesh->GetNumberOfPoints();
    }
    else {
        number_of_samples = res[0]*res[1]*res[2];
    }
    float *data;
    int nscalars = 0;
    int offset = 0;
    if (v!=nullptr) {
        nscalars += 3;
    }
    if (s != nullptr) {
        nscalars += 1;
        offset = 1;
    }
    data = (float*)calloc(number_of_samples*nscalars, sizeof(float));

    timer.start();
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (size_t n=0 ; n<number_of_samples ; ++n) {
            spurt::fvec3 p;
            if (has_geometry) {
                double pp[3];
                target_mesh->GetPoint(n, pp);
                std::copy(pp, pp+3, p.begin());
            }
            else {
                int i = n % res[0];
                int j = n / res[0];
                int k = j / res[1];
                j = j % res[1];
                p = bounds.min() + step*spurt::fvec3(i, j, k);
            }
            if (vintp != nullptr) (*vintp)(0, &p[0], &data[nscalars*n+offset]);
            if (sintp != nullptr) (*sintp)(0, &p[0], &data[nscalars*n]);

            int thread = 0;
#if _OPENMP
            thread = omp_get_thread_num();
#endif
            if (!thread) {
                double elapsed = _timer.elapsed();
                std::cout << "\r" << n << " / " << number_of_samples
                          << " (" << 100.*(double)n/(double)number_of_samples << "\%) in "
                          << elapsed << " s. (" << (double)n/elapsed << " Hz)         \r"
                          << std::flush;
            }
        }
    }
    std::cout << '\n';

    if (has_geometry) {
        if (s != nullptr) {
            float* scalars = new float[number_of_samples];
            for (size_t i=0; i<number_of_samples; ++i) {
                scalars[i] = data[nscalars*i];
            }
            vtk_utils::add_scalars_from_carray(target_mesh, scalars, number_of_samples, true, "pressure", true);
            delete[] scalars;
        }
        if (v != nullptr) {
            float* vectors = new float[3*number_of_samples];
            for (size_t i=0; i<number_of_samples; ++i) {
                vectors[3*i  ] = data[nscalars*i+offset  ];
                vectors[3*i+1] = data[nscalars*i+offset+1];
                vectors[3*i+2] = data[nscalars*i+offset+2];
            }
            vtk_utils::add_vectors_from_carray<float, VTK_SMART(vtkDataSet), 3>(target_mesh, vectors, number_of_samples, true, "velocity", true);
            delete[] vectors;
        }
        vtk_utils::saveVTK(target_mesh, base(out_name) + '.' + extension(grid_name));
    }
    else {
        // NRRD file storage
        spurt::nrrd_utils::nrrd_params<float, 4> params;
        params.mins()[0] = 0;
        params.mins()[1] = bounds.min()[0];
        params.mins()[2] = bounds.min()[1];
        params.mins()[3] = bounds.min()[2];
        params.spacings()[0] = 1;
        params.spacings()[1] = step[0];
        params.spacings()[2] = step[1];
        params.spacings()[3] = step[2];
        params.sizes()[0] = nscalars;
        params.sizes()[1] = res[0];
        params.sizes()[2] = res[1];
        params.sizes()[3] = res[2];
        params.labels()[0] = std::string((sintp != nullptr) ? "pressure;" : "") + std::string((vintp != nullptr) ? "Vx;Vy;Vz" : "");
        params.labels()[1] = "x";
        params.labels()[2] = "y";
        params.labels()[3] = "z";

        spurt::nrrd_utils::writeNrrdFromParams(data, out_name, params);
    }
    delete[] data;
    delete m;
    delete s;
    delete v;
    std::cout << out_name << " exported\n";
    return 0;
}
