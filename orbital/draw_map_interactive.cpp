#include <iostream>
#include <fstream>
#include <vector>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>

#include "new_universal_poincare_map.hpp"
#include "cr3bp.hpp"

#include <VTK/vtk_utils.hpp>
#include <misc/sort.hpp>
#include <graphics/colors.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

double K;
nvis::vec2 res;
int maxi, nsamples;
std::string filename, input_name;
double eps;
double C;
double mu;
int seeding;
nvis::bbox2 bounds;
bool verbose = false;
bool show_boundary=false;
bool show_edges=false;
int connection_method=0;
int scale_order=3;

typedef nvis::fixed_vector< double, 6 > vec6;

double mean_dist(const std::vector<nvis::vec2>& x, int p) {
    double d=0;
    for (int i=0; i<x.size()-p; ++i) {
        d+=nvis::norm(x[i+p]-x[i]);
    }
    return d/(x.size()-p);
}

int best_period(double& _d, const std::vector<nvis::vec2>& x, int maxp) {
    int bestp=-1;
    _d=std::numeric_limits<double>::max();
    // std::cout << "maxp=" << maxp << '\n';
    for (int p=1; p<=maxp ; ++p) {
        // std::cout << "x.size()=" << x.size() << '\n';
        double d = mean_dist(x, p);
        if (d<_d) {
            bestp=p;
            _d=d;
        }
    }
    return bestp;
}

double curvature(const std::vector<nvis::vec2>& pts, int period) {
    int n=0;
    double total_curvature=0;
    for (int i=0; i<period; ++i) {
        int j;
        for (j=i; j+2*period<pts.size(); j+=period, ++n) {
            const nvis::vec2& p0=pts[j];
            const nvis::vec2& p1=pts[j+period];
            const nvis::vec2& p2=pts[j+2*period];
            total_curvature += nvis::inner(p1-p0, p2-p1)/nvis::norm(p1-p0)/nvis::norm(p2-p1);
        }
    }
    return total_curvature/(double)n;
}

void generate_seeds(std::vector< nvis::vec2 >& p)
{
    const double& minx = bounds.min()[0];
    const double& miny = bounds.min()[1];
    const double& maxx = bounds.max()[0];
    const double& maxy = bounds.max()[1];

    p.clear();
    switch (seeding) {
    case 0: {
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(minx + drand48()*(maxx - minx),
                                   miny + drand48()*(maxy - miny)));
        break;
    }
    case 1: {
        double dx = (maxx - minx) / (double)nsamples;
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(minx + (double)n*dx, 0.5*(miny + maxy)));
        break;
    }
    case 2: {
        double dy = (maxy - miny) / (double)nsamples;
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(0.5*(minx + maxx), miny + (double)n*dy));
        break;
    }
    case 3: {
        double dx = (maxx - minx) / (double)nsamples;
        double dy = (maxy - miny) / (double)nsamples;
        for (int n = 0 ; n < nsamples ; ++n)
            p.push_back(nvis::vec2(minx + (double)n*dx, miny + (double)n*dy));
        break;
    }
    default: {
        std::cout << "unknown seeding pattern. using random sampling" << std::endl;
        seeding = 0;
        generate_seeds(p);
    }
    }
}

std::string me;
void displayUsageAndExit(const std::string& what) {
    if (what.size()) std::cerr << "ERROR: " << what << '\n';
    std::cerr
        << "USAGE: " << me << " [options]\n"
        << "DESCRIPTION: compute Poincare plot of circular restricted 3-body problem\n"
        << "OPTIONS:\n"
        << "   -h | --help                  Print this information\n"
        << "   -o | --output <string>       Output file name\n"
        << "   -r | --res <int> (x2)        Resolution\n"
        << "   -i | --iterations <int>      Number of iterations\n"
        << "   -b | --bounds <float> (x4)   Computation bounds\n"
        << "   -e | --eps <float>           Integration precision\n"
        << "   -n | --samples <int>         Number of sample trajectories\n"
        << "      | --system <string>       System to consider\n"
        << "   -C <float>                   C constant\n"
        << "   -m | --mu <float>            mu constant\n"
        << "   -s | --seed <int>            Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal)\n"
        << "      | --show_boundary         Indicate boundary of selected region by a quad\n"
        << "      | --show_edges            Connect points of Poincare plot\n"
        << "      | --connection            Method to determine connectivity\n"
        << "      | --scale_order           Order of color scale progression\n"
        << "   -v | --verbose               Select verbose output\n";
    exit(1);
}

void readParams(const std::string& filename) {
    std::fstream input(filename, std::ios::in);
    while (!input.eof()) {
        std::string tmp, key;
        char c;
        input >> key;
        if (key == "mu") input >> mu;
        else if (key == "C") input >> C;
        else if (key == "xmin") input >> bounds.min()[0];
        else if (key == "ymin") input >> bounds.min()[1];
        else if (key == "xmax") input >> bounds.max()[0];
        else if (key == "ymax") input >> bounds.max()[1];
        else if (key == "eps") input >> eps;
        else if (key == "iter") input >> maxi;
        else if (key == "number") input >> nsamples;
        else if (key == "res") input >> res[0] >> res[1];
    }
}

void sort(double& d, std::vector<nvis::ivec2>& edges, const std::vector<nvis::vec2>& pts) {
    const int N=pts.size();
    edges.clear();
    d=0;
    for (int i=0; i<N; ++i) {
        std::vector<double> dist(N);
        for (int j=0; j<N; ++j) {
            dist[j] = nvis::norm(pts[i]-pts[j]);
        }
        std::vector<unsigned int> sorted;
        xavier::sort(dist, sorted);
        edges.push_back(nvis::ivec2(i, sorted[1]));
        edges.push_back(nvis::ivec2(i, sorted[2]));
        d+=nvis::norm(pts[i]-pts[sorted[1]]);
        d+=nvis::norm(pts[i]-pts[sorted[2]]);
    }
    d /= (2*N);
}

struct Less {
    template<typename T>
    bool operator()(const std::pair<double, T>& p1, const std::pair<double, T>& p2) {
        return p1.first < p2.first;
    }
};

int main(int argc, char* argv[])
{
    me = argv[0];
    bounds.min()[0] = -0.6985;
    bounds.min()[1] = -1.3335;
    bounds.max()[0] = -0.5080;
    bounds.max()[1] = 1.3335;
    C = 3.000;
    mu = Earth_Moon_mu;
    res = nvis::vec2(512, 512);
    eps = 1.0e-7;
    nsamples = 100;
    maxi = 1000;
    filename = "none";
    input_name = "none";

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") displayUsageAndExit("");
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) displayUsageAndExit("missing output name");
            filename = argv[++i];
        }
        else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) displayUsageAndExit("missing resolution");
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
        }
        else if (arg == "-C") {
            if (i == argc-1) displayUsageAndExit("missing C constant");
            C = atof(argv[++i]);
        }
        else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) displayUsageAndExit("missing mu constant");
            mu = atof(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) displayUsageAndExit("missing bounds");
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
        }
        else if (arg == "-i" || arg == "--iterations") {
            if (i == argc-1) displayUsageAndExit("missing number of iterations");
            maxi = atof(argv[++i]);
        }
        else if (arg == "-n" || arg == "--samples") {
            if (i == argc-1) displayUsageAndExit("missing number of samples");
            nsamples = atoi(argv[++i]);
        }
        else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) displayUsageAndExit("missing precision");
            eps = atof(argv[++i]);
        }
        else if (arg == "-s" || arg == "--seeding") {
            if (i == argc-1) displayUsageAndExit("missing precision");
            seeding = atoi(argv[++i]);
        }
        else if (arg == "-f" || arg == "--file") {
            if (i == argc-1) displayUsageAndExit("missing input filename");
            readParams(argv[++i]);
        }
        else if (arg == "--system") {
            if (i == argc-1) displayUsageAndExit("missing system name");
            std::string s(argv[++i]);
            if (s=="earth-moon" || s=="Earth-Moon") {
                mu = Earth_Moon_mu;
            }
            else if (s=="jupiter-europa" || s=="Jupiter-Europa") {
                mu = Jupiter_Europa_mu;
            }
        }
        else if (arg == "--show_boundary") {
            show_boundary=true;
        }
        else if (arg == "--show_edges") {
            show_edges=true;
        }
        else if (arg == "-v" || arg == "--verbose") {
            verbose=true;
        }
        else if (arg == "--connection") {
            if (i == argc-1) displayUsageAndExit("missing connectivity name");
            std::string s(argv[++i]);
            if (s=="period" || s=="periodic") {
                connection_method=0;
            }
            else if (s=="distance" || s=="dist" || s=="mindist" || s=="min_dist" || s=="min_distance") {
                connection_method=1;
            }
            else {
                std::cout << "Warning: connectivity name (" << s << ") unrecognized\n";
            }
        }
        else if (arg == "--scale_order") {
            if (i == argc-1) displayUsageAndExit("missing scale_order");
            scale_order = atoi(argv[++i]);
        }
        else {
            displayUsageAndExit("unrecognized argument");
        }
    }

    if (nvis::all(bounds.min() == bounds.max()))
        displayUsageAndExit("Missing bounds");

    if (verbose) {
        std::cout << "parameters: resolution = " << res << '\n';
        std::cout << "bounds=\n" << bounds << std::endl;
    }

    if (nsamples == -1) nsamples = res[0];
    if (maxi == -1) maxi = res[0];

    float *hits = (float*)calloc(3 * res[0] * res[1], sizeof(float));

    srand48(987654321);

    VTK_CREATE(vtkRenderer, renderer);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);

    std::vector< nvis::vec2 > seeds;
    generate_seeds(seeds);

    unsigned int count = 0;
    unsigned int pct, last_pct = 0;

    cr3bp dummy;
    planar_section<6>* section = new planar_section<6>(dummy.plane());
    int failed = 0;

    std::vector<nvis::fvec3> scale;
    xavier::spiral_scale(scale, 10, 0.1, 1, 1, 0, scale_order);

    std::vector<std::pair<double, VTK_SMART(vtkPolyData)> > polydatas;
    std::vector<double> values;

    std::fstream csv("tmp.csv", std::ios::out);

    std::vector< std::vector<double> > all_times;

    nvis::timer _timer;
#pragma omp parallel
    {
#ifdef _OPENMP
        const int thread = omp_get_thread_num();
#endif
        cr3bp *field = new cr3bp(C, mu);
        new_universal_poincare_map<cr3bp, 6, planar_section<6> > _map(field, section);
        _map.precision(eps);

#pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < nsamples ; ++n) {
            ++count;
            nvis::vec2 x = seeds[n];
            vec6 y;
            try {
                y = field->unproject(x);
            }
            catch (...) {
                if (verbose) std::cout << "could not start at " << x << std::endl;
                ++failed;
                continue;
            }

            if (verbose) std::cout << "starting at " << x << '\n';

            std::vector< vec6 > __hits;
            std::vector< double > __times;
            try {
                _map.map(y, __hits, __times, maxi);
                all_times.push_back(__times);
            }
            catch (...) {
                std::cout << "caught exception at " << x << std::endl;
            }

            std::vector<nvis::vec2> points;
            points.push_back(x);
            for (int i=0; i<__hits.size(); ++i) {
                nvis::vec2 p = field->project(__hits[i]);
                points.push_back(p);
            }

            VTK_SMART(vtkPolyData) poly=vtk_utils::make_points(points);
            vtk_utils::add_vertices(poly);

            double meandist;
            if (true) {
                VTK_CREATE(vtkCellArray, cells);
                if (connection_method == 0) {
                    int period=best_period(meandist, points, points.size()/10);
                    for (int i=0; i<points.size()-points.size()/10; ++i) {
                        cells->InsertNextCell(2);
                        cells->InsertCellPoint(i);
                        cells->InsertCellPoint(i+period);
                    }
                    csv << meandist << ", " << curvature(points, period)
                        << ", " << x[0] << ", " << x[1] << '\n';
                    std::cout << meandist << " / "
                        << curvature(points, period) << '\n';
                } else {
                    std::vector<nvis::ivec2> edges;
                    sort(meandist, edges, points);
                    for (int i=0; i<edges.size(); ++i) {
                        cells->InsertNextCell(2);
                        cells->InsertCellPoint(edges[i][0]);
                        cells->InsertCellPoint(edges[i][1]);
                    }
                }
                if (show_edges) poly->SetLines(cells);
            }
            values.push_back(meandist);

            polydatas.push_back(std::pair<double, VTK_SMART(vtkPolyData)>(meandist, poly));

            pct = floor((double)count / (double)nsamples * 100.);
#ifdef _OPENMP
            if (thread == 0 && pct > last_pct)
#else
            if (pct > last_pct)
#endif
            {
                std::ostringstream os;
                os << count << " curves computed so far (" << pct << "%)\n";
                std::cout << os.str();
                last_pct = pct;
            }
        }
    }
    double delta_t = _timer.elapsed();
    std::cout << '\n';
    csv.close();

    std::cout << failed << " integrations failed\n";
    std::cout << "integration of " << nsamples << " orbits took " << delta_t << " s. ("
         << static_cast<double>(nsamples)/delta_t << " Hz)\n";

    std::cout << "\nTimes associated with each iterations:\n";
    for (auto it=all_times.begin(); it!=all_times.end(); ++it) {
        auto times = *it;
        std::copy(times.begin(), times.end(), std::ostream_iterator<double>(std::cout, ", "));
        std::cout << '\n';
    }

    xavier::adaptive_color_map<double> __cmap(values, scale, false);
    VTK_CREATE(vtkColorTransferFunction, cmap);
    for (int i=0; i<__cmap.t.size(); ++i) {
        nvis::fvec3 c=__cmap.t[i];
        cmap->AddRGBPoint(__cmap.t[i], c[0], c[1], c[2]);
    }

    std::sort(polydatas.begin(), polydatas.end(), Less());
    double minv = polydatas.front().first;
    double maxv = polydatas.back().first;

    for (int i=0; i<polydatas.size(); ++i) {
        double v = polydatas[i].first;
        double t = (v-minv)/(maxv-minv);
        double u = pow(t, scale_order);
        nvis::fvec3 col=__cmap(polydatas[i].first);
        VTK_MAKE_ACTOR(actor, polydatas[i].second);
        actor->GetProperty()->SetColor(col[0], col[1], col[2]);
        actor->GetProperty()->SetPointSize(2*t);
        renderer->AddActor(actor);
    }
    // for (int i=0; i<polydatas.size(); ++i) {
    //     int N=polydatas[i].second->GetNumberOfPoints();
    //     std::vector<double> values(N, polydatas[i].first);
    //     vtk_utils::add_scalars(polydatas[i].second, values);
    //     VTK_SMART(vtkPolyData) pd = vtk_utils::make_spheres(polydatas[i].second,
    //         *std::max_element(values.begin(), values.end()), 12, 12, true);
    //     VTK_MAKE_ACTOR(actor, pd);
    //     actor->GetMapper()->ScalarVisibilityOn();
    //     actor->GetMapper()->SetLookupTable(cmap);
    //     renderer->AddActor(actor);
    // }

    if (show_boundary) {
        std::vector<nvis::vec2> corners(4);
        corners[0] = bounds.min();
        corners[2] = bounds.max();
        corners[1][0] = corners[2][0];
        corners[1][1] = corners[0][1];
        corners[3][0] = corners[0][0];
        corners[3][1] = corners[2][1];
        VTK_SMART(vtkPolyData) frame=vtk_utils::make_points(corners);

        VTK_CREATE(vtkCellArray, quad);
        quad->InsertNextCell(4);
        quad->InsertCellPoint(0);
        quad->InsertCellPoint(1);
        quad->InsertCellPoint(2);
        quad->InsertCellPoint(3);
        frame->SetLines(quad);

        VTK_MAKE_ACTOR(frame_actor, frame);
        frame_actor->GetProperty()->SetEdgeVisibility(1);
        frame_actor->GetProperty()->SetLineWidth(1);
        frame_actor->GetProperty()->SetColor(1,1,1);
        renderer->AddActor(frame_actor);
    }

    interactor->Initialize();
    renderer->GetActiveCamera()->SetParallelProjection(1);
    nvis::vec2 center=bounds.center();
    renderer->GetActiveCamera()->SetPosition(center[0], center[1], 1);
    renderer->GetActiveCamera()->SetFocalPoint(center[0], center[1], 0);
    renderer->GetActiveCamera()->SetViewUp(0, 1, 0);
    renderer->GetActiveCamera()->SetParallelScale(0.5*(bounds.max()[1]-bounds.min()[1]));
    window->Render();
    interactor->Start();

    return 0;
}
