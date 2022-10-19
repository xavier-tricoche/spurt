#include <iostream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <Eigen/SVD>
#include <math/MLS.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <sfcnn.hpp>

#include <teem/nrrd.h>

#include <vtkGenericCell.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>
#include <vtkDataArray.h>

#include <vtkCellTree.hpp>
#include <vtkDLRReader.hpp>
#include <vtk_mesh_traits.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
int fit_order;
int res[3];
int n_order;

typedef sfcnn<nvis::fvec3, 3, float> tree_type;

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
    hestOptAdd(&hopt, "i",      "input",                airTypeString,  1, 1, &name_in,             NULL,       "input file name");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &name_out,            NULL,       "output file name (NRRD)");
    hestOptAdd(&hopt, "r",      "resolution",           airTypeInt,     3, 3, &res,                 NULL,       "resampling resolution");
    hestOptAdd(&hopt, "k",      "order",                airTypeInt,     0, 1, &fit_order,           "0",        "resampling order");
    hestOptAdd(&hopt, "n",      "# neighbors",          airTypeInt,     0, 1, &n_order,             "4",        "number of neighbors involved in fit");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Resample input unstructured mesh over regular grid using smooth reconstruction kernels",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

vtkDataSet* read(const std::string& filename)
{

    std::string prefix = filename.substr(0, filename.rfind('.'));
    std::string suffix = filename.substr(filename.rfind('.'), filename.size()-prefix.size()-1);
    
    std::cerr << "input file prefix: " << prefix << std::endl;
    std::cerr << "input file suffix: " << suffix << std::endl;
    
    vtkDataSet* data;
    if (!strcmp(suffix.c_str(), "vtk")) {
        vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
        reader->SetFileName(name_in);
        reader->Update();
        data = reader->GetOutput();
        data->Register(0);
        reader->Delete();
    } else if (!strcmp(suffix.c_str(), "grid") || !strcmp(suffix.c_str(), "pval")) {
        vtkDLRReader* reader = vtkDLRReader::New();
        reader->SetFileName(name_in);
        reader->Update();
        data = reader->GetOutput();
        data->Register(0);
        reader->Delete();
    } else {
        std::cerr << "unrecognized file extension: " << suffix << std::endl;
        exit(-1);
    }
    
    return data;
}

nvis::fvec3* all_pts;
tree_type* create_tree(vtkDataSet* data)
{
    int npts = data->GetNumberOfPoints();
    all_pts = new nvis::fvec3[npts];
    nvis::vec3 x;
    for (int i=0 ; i<npts ; ++i) {
        data->GetPoint(i, x.begin());
        all_pts[i] = x;
    }
    
    tree_type* tree = new tree_type(all_pts, 3);
    
    return tree;
}

static int nbchannels[3][3] = { { 1, 1, 1 }, { 1, 2, 4 }, { 1, 3, 9 } };

template<typename value_type>
void fit_one(Eigen::MatrixXf& res, unsigned int id,
             vtkDataSet* data, vtkDataArray* vals,
             tree_type* tree)
{

    MLS::weighted_least_squares<value_type, nvis::vec3> WLS(3, fit_order, value_type::size());
    
    nvis::vec3 xd, y;
    data->GetPoint(id, &(xd[0]));
    nvis::fvec3 xf = xd;
    std::vector<long unsigned int> neighbors;
    std::vector<double> distances;
    tree->ksearch(xf, n_order+1, neighbors, distances);
    double radius = *std::max_element(distances.begin(), distances.end());
    
    std::vector<nvis::vec3> points(neighbors.size());
    std::vector<value_type> attributes(neighbors.size());
    for (int n=0 ; n<neighbors.size() ; ++n) {
        data->GetPoint(neighbors[n], &(points[n][0]));
        vals->GetTuple(neighbors[n], &(attributes[n][0]));
    }
    
    Eigen::MatrixXd r;
    WLS(r, points, attributes, xd, radius);
    res = r.base().cast<float>();
}

inline void fit_one(Eigen::MatrixXf& res, unsigned int id,
                    vtkDataSet* data, vtkDataArray* vals,
                    tree_type* tree, int order)
{

    int size = nbchannels[2][order];
    switch (size) {
        case 1:
            return fit_one<nvis::fixed_vector<double, 1> >(res, id, data, vals, tree);
        case 2:
            return fit_one<nvis::fixed_vector<double, 2> >(res, id, data, vals, tree);
        case 3:
            return fit_one<nvis::fixed_vector<double, 3> >(res, id, data, vals, tree);
        case 4:
            return fit_one<nvis::fixed_vector<double, 4> >(res, id, data, vals, tree);
        case 9:
            return fit_one<nvis::fixed_vector<double, 9> >(res, id, data, vals, tree);
    }
}

void fit(std::vector<Eigen::MatrixXf>& values, vtkDataSet* data, tree_type* tree)
{

    size_t npts = data->GetNumberOfPoints();
    vtkDataArray* rhs;
    int order;
    if (data->GetPointData()->GetScalars()) {
        order = 0;
        rhs = data->GetPointData()->GetScalars();
    } else if (data->GetPointData()->GetVectors()) {
        order = 1;
        rhs = data->GetPointData()->GetVectors();
    } else if (data->GetPointData()->GetTensors()) {
        order = 2;
        rhs = data->GetPointData()->GetTensors();
    } else {
        std::cerr << "unrecognized data attribute\n";
        exit(-1);
    }
    
    int dof = MLS::dof(3, fit_order);
    int nbc = nbchannels[3][order];
    values.resize(npts);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int i=0 ; i<npts ; ++i) {
            fit_one(values[i], i, data, rhs, tree, order);
        }
    }
}

void resample(const nvis::bbox3& bounds,
              const std::vector<Eigen::MatrixXf>& values,
              vtkCellTree* locator)
{
    int size = values[0].size();
    nvis::ivec3 dims(res[0], res[1], res[2]);
    xavier::RasterGrid<3> grid(dims, bounds);
    int npts = grid.size();
    float* data = (float*)calloc(npts*size, sizeof(float));
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int i=0 ; i<npts ; ++i) {
            nvis::vec3 x = grid(grid.coord(i));
            double tol2, pcoords[3], weights[8];
            vtkGenericCell* cell;
            vtkIdType cellid = locator->FindCell(&(x[0]), tol2, cell, pcoords, weights);
            if (cellid != -1) {
                Eigen::MatrixXf m = weights[0] * values[cell->GetPointId(0)];
                for( int j=1; j<cell->GetNumberOfPoints(); ++j ) {
                    m += weights[j]*values[cell->GetPointId(j)];
                }
                
                int id = i*size;
                for (int r=0 ; r<m.rows() ; ++r)
                    for (int c=0 ; c<m.cols() ; ++c) {
                        data[id++] = m(r, c);
                    }
            }
        }
    }
    
    std::vector<double> spc(4);
    std::vector<size_t> sz(4);
    sz[0] = size;
    spc[0] = airNaN();
    for (int i=0 ; i<3 ; ++i) {
        sz[i+1] = res[i];
        spc[i+1] = grid.step()[i];
    }
    xavier::writeNrrd(data, name_out, nrrdTypeFloat, sz, spc);
}


int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    typedef mesh_traits<vtkDataSet*>    mtraits;
    
    std::string filename(name_in);
    vtkDataSet* data = read(filename);
    std::cerr << filename << " successfully imported\n";
    
    tree_type* tree = create_tree(data);
    std::cerr << "kD tree constructed\n";
    
    std::vector<Eigen::MatrixXf> values;
    fit(values, data, tree);
    std::cerr << "polynomial fits computed\n";
    
    vtkCellTree* locator = vtkCellTree::New();
    locator->SetDataSet(data);
    locator->LazyEvaluationOff();
    locator->CacheCellBoundsOff();
    locator->BuildLocator();
    
    double b[6];
    data->GetBounds(b);
    nvis::bbox3 bounds;
    bounds.min() = nvis::vec3(b[0], b[2], b[4]);
    bounds.max() = nvis::vec3(b[1], b[3], b[5]);
    resample(bounds, values, locator);
    
    exit(0);
}
