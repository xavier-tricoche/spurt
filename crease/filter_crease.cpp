#include <vtk/vtk_utils.hpp>
#include <vtk/filter_ccs.hpp>
#include <misc/cxxopts.hpp>
#include <string>
#include <math/types.hpp>
#include "vtkTransformPolyDataFilter.h"
#include <math/stat.hpp>
using namespace spurt;

constexpr double minus_infinity = std::numeric_limits<double>::lowest();
constexpr double plus_infinity = std::numeric_limits<double>::max();
constexpr double epsilon = std::numeric_limits<double>::min();
constexpr int very_large = std::numeric_limits<int>::max();

VTK_SMART(vtkPolyData) translate(VTK_SMART(vtkPolyData) input, const vec3& t)
{
    if (spurt::norm(t) == 0) return input;
    VTK_CREATE(vtkTransform, translate);
    translate->Identity();
    translate->Translate(t[0], t[1], t[2]);
    VTK_CREATE(vtkTransformPolyDataFilter, filter);
    filter->SetTransform(translate);
    filter->SetInputData(input);
    filter->Update();
    return filter->GetOutput();
}

std::vector<double> array2vector(VTK_SMART(vtkAbstractArray) _array) {
    VTK_SMART(vtkDataArray) array = vtkDataArray::SafeDownCast(_array);
    std::vector<double> values;
    long ntuples = array->GetNumberOfTuples();
    int ncomps = array->GetNumberOfComponents();
    std::cout << "there are " << ntuples << " tuples and " << ncomps << " comps\n";
    for (long i=0; i<ntuples; ++i) {
        double* ptr = array->GetTuple(i);
        if (ncomps == 1) {
            values.push_back(*ptr);
        }
        else {
            values.push_back(sqrt(std::inner_product(ptr, ptr+ncomps, ptr, 0)));
        }
    }
    std::cout << "values contains " << values.size() << " entries\n";
    return values;
}

void array2vector2(std::vector<double>& values, VTK_SMART(vtkAbstractArray) _array) {
    VTK_SMART(vtkDataArray) array = vtkDataArray::SafeDownCast(_array);
    long ntuples = array->GetNumberOfTuples();
    int ncomps = array->GetNumberOfComponents();
    std::cout << "there are " << ntuples << " tuples and " << ncomps << " comps\n";
    values.resize(ntuples);
    for (long i=0; i<ntuples; ++i) {
        double* ptr = array->GetTuple(i);
        if (ncomps == 1) {
            values[i] = (*ptr);
        }
        else {
            values[i] = sqrt(std::inner_product(ptr, ptr+ncomps, ptr, 0));
        }
    }
    std::cout << "values contains " << values.size() << " entries\n";
}

VTK_SMART(vtkPolyData) filter_cells(VTK_SMART(vtkPolyData) input, 
                                    const std::vector<bool>& selected)
{
    std::cout << "entering filter_cells\n";
    VTK_SMART(vtkCellArray) vertices = input->GetVerts();
    VTK_SMART(vtkCellArray) lines = input->GetLines();
    VTK_SMART(vtkCellArray) polys = input->GetPolys();
    int nverts = input->GetNumberOfVerts();
    int nlines = input->GetNumberOfLines();
    int npolys = input->GetNumberOfPolys();
    assert(nverts + nlines + npolys == selected.size());

    VTK_CREATE(vtkCellArray, newverts);
    VTK_CREATE(vtkCellArray, newlines);
    VTK_CREATE(vtkCellArray, newpolys);

    VTK_SMART(vtkCellData) celldata = input->GetCellData();
    int natts = celldata->GetNumberOfArrays();
    std::cout << "there are " << natts << " cell attributes\n";
    std::vector< VTK_SMART(vtkDataArray) > out_attributes(natts);
    for (int i=0; i<natts; ++i) {
        VTK_SMART(vtkDataArray) data = celldata->GetArray(i);
        std::cout << "attribute #" << i << " is called " << data->GetName() << '\n';
        std::cout << "this data set contains " << data->GetNumberOfTuples() << " tuples and each tuple contains " << data->GetNumberOfComponents() << " coefficients\n";
        out_attributes[i] = vtkDataArray::CreateDataArray(data->GetDataType());
        out_attributes[i]->SetName(data->GetName());
        out_attributes[i]->SetNumberOfComponents(data->GetNumberOfComponents());
        std::cout << "copy attribute " << i << " is called " << out_attributes[i]->GetName() << '\n';
        std::cout << "It contains " << out_attributes[i]->GetNumberOfTuples() << " tuples and each tuple contains " << out_attributes[i]->GetNumberOfComponents() << " coefficients\n";
    }

    // Note: what about strips?
    for (int what=0; what<3; ++what)
    {
        VTK_SMART(vtkCellArray) old_cells;
        VTK_SMART(vtkCellArray) new_cells;
        int offset = 0;
        if (what == 0) { // Verts
            old_cells = vertices;
            new_cells = newverts;
        }
        else if (what == 1) { // Lines
            old_cells = lines;
            offset = nverts;
            new_cells = newlines;
        } 
        else { // Polys
            old_cells = polys;
            offset = nverts + nlines;
            new_cells = newpolys;
        }
        VTK_CREATE(vtkIdList, alist);
        int cellid = offset;
        for (auto it=old_cells->NewIterator(); 
             !it->IsDoneWithTraversal(); it->GoToNextCell(), ++cellid) {
            if (selected[cellid]) {
                it->GetCurrentCell(alist);
                new_cells->InsertNextCell(alist);
                for (int i=0; i<natts; ++i) {
                    VTK_SMART(vtkDataArray) data = celldata->GetArray(i);
                    out_attributes[i]->InsertNextTuple(data->GetTuple(cellid));
                }
            }
        }
    }

    VTK_CREATE(vtkPolyData, output);
    output->DeepCopy(input);
    output->SetPolys(newpolys);
    output->SetLines(newlines);
    output->SetVerts(newverts);
    for (int i=0; i<natts; ++i) {
        output->GetCellData()->RemoveArray(i);
    }
    for (int i=0; i<natts; ++i) {
        output->GetCellData()->AddArray(out_attributes[i]);
    }
    std::cout << "leaving filter_cells: output contains " 
              << output->GetNumberOfPoints() << " points and " 
              << output->GetNumberOfCells() << " cells\n";
    
    return output;
}

VTK_SMART(vtkPolyData) filter_points(VTK_SMART(vtkPolyData) input, 
                                     const std::vector<bool>& selected)
{
    std::cout << "entering filter_points\n";
    assert(selected.size() == input->GetNumberOfPoints());
    std::vector<bool> selected_cells(input->GetNumberOfCells(), false);
    VTK_CREATE(vtkIdList, alist);
    int nselected = 0;
    for (int cellid=0; cellid<selected_cells.size(); ++cellid) 
    {
        input->GetCellPoints(cellid, alist);
        bool included = true;
        for (const vtkIdType& id : *alist) {
            if (!selected[id]) {
                included = false;
                break;
            }
        }
        selected_cells[cellid] = included;
        if (selected_cells[cellid]) ++nselected;
    }
    std::cout << "leaving filter_points\n";
    std::cout << "There were " << nselected << " selections\n";
    return filter_cells(input, selected_cells);
}

template<typename T>
VTK_SMART(vtkPolyData) filter_by_value(VTK_SMART(vtkPolyData) input,
                                       const std::string& array_name, 
                                       T min, T max, bool cellwise)
{
    std::cout << "\n\nFiltering a dataset with " 
        << input->GetNumberOfPoints() << " points and " 
        << input->GetNumberOfCells() << " cells by " 
        << array_name << " " << (cellwise ? "cellwise" : "pointwise") << "\n\n";
    std::cout << "value range is " << min << " -> " << max << '\n';

    VTK_SMART(vtkDataArray) values;
    std::cout << "getting data array... " << std::flush;
    if (cellwise) {
        values = input->GetCellData()->GetArray(array_name.c_str());
    }
    else {
        values = input->GetPointData()->GetArray(array_name.c_str());
    }
    std::cout << "done\n";
    std::vector<bool> selected(values->GetNumberOfTuples(), false);
    std::cout << "checking individual tuples... " << std::flush;
    int nselected = 0;
    for (int i=0; i<selected.size(); ++i)
    {
        T v = values->GetTuple1(i);
        selected[i] = (v >= min) && (v <= max);
        if (selected[i]) ++nselected;
    }
    std::cout << "done\n";
    std::cout << "There were " << nselected << " selections\n";
    if (cellwise) return filter_cells(input, selected);
    else return filter_points(input, selected);
}

VTK_SMART(vtkPolyData) shrink(VTK_SMART(vtkPolyData) input)
{
    std::map<long, long> old2new;
    std::vector<long> new2old;

    VTK_SMART(vtkCellArray) vertices = input->GetVerts();
    VTK_SMART(vtkCellArray) lines = input->GetLines();
    VTK_SMART(vtkCellArray) polys = input->GetPolys();
    VTK_CREATE(vtkCellArray, newverts);
    VTK_CREATE(vtkCellArray, newlines);
    VTK_CREATE(vtkCellArray, newpolys);
    VTK_CREATE(vtkPolyData, shrunk);

    // loop over all cells of all cell kinds and reindex vertices
    for (int what=0; what<3; ++what)
    {
        VTK_SMART(vtkCellArray) old_cells;
        VTK_SMART(vtkCellArray) new_cells;
        if (what == 0) {
            old_cells = vertices;
            new_cells = newverts;
        }
        else if (what == 1) {
            old_cells = lines;
            new_cells = newlines;
        }
        else if (what == 2) {
            old_cells = polys;
            new_cells = newpolys;
        }

        VTK_CREATE(vtkIdList, alist);
        for (auto it=old_cells->NewIterator(); 
             !it->IsDoneWithTraversal(); it->GoToNextCell()) {
            it->GetCurrentCell(alist);
            long nids = alist->GetNumberOfIds();
            for (long n=0; n<nids; ++n) {
                long id = alist->GetId(n);
                auto iter = old2new.find(id);
                if (iter == old2new.end()) {
                    old2new[id] = new2old.size();
                    alist->SetId(n, new2old.size());
                    new2old.push_back(id);
                }
                else {
                    alist->SetId(n, iter->second);
                }
            }
            new_cells->InsertNextCell(alist);
        }
    }
    shrunk->SetVerts(newverts);
    shrunk->SetLines(newlines);
    shrunk->SetPolys(newpolys);

    // Update points array
    VTK_SMART(vtkDataArray) oldpos = input->GetPoints()->GetData();
    VTK_CREATE(vtkFloatArray, newpos);
    newpos->SetNumberOfComponents(3);
    newpos->SetNumberOfTuples(new2old.size());
    for (long i=0; i<new2old.size(); ++i) {
        auto ptr = oldpos->GetTuple3(new2old[i]);
        newpos->SetTuple3(i, ptr[0], ptr[1], ptr[2]);
    }
    VTK_CREATE(vtkPoints, newpoints);
    newpoints->SetData(newpos);
    shrunk->SetPoints(newpoints);
    auto nnew = shrunk->GetNumberOfPoints();
    auto nold = input->GetNumberOfPoints();
    std::cout << "shrunk dataset contains " << nnew << " points (from " << nold << ")\n";
    std::cout << "ratio: " << (double)nnew/(double)nold*100. << "%\n";

    // Update point attributes
    VTK_SMART(vtkPointData) pointdata = input->GetPointData();
    int natts = pointdata->GetNumberOfArrays();
    std::cout << "there are " << natts << " points attributes\n";
    VTK_SMART(vtkPointData) newpointdata = shrunk->GetPointData();
    for (int n=0; n<natts; ++n) {
        VTK_SMART(vtkDataArray) data = pointdata->GetArray(n);
        std::cout << "attribute #" << n << " is called " << data->GetName() << '\n';
        std::cout << "this data set contains " << data->GetNumberOfTuples() << " tuples and each tuple contains " << data->GetNumberOfComponents() << " coefficients\n";
        VTK_SMART(vtkDataArray) att = vtkDataArray::CreateDataArray(data->GetDataType());
        att->SetName(data->GetName());
        att->SetNumberOfComponents(data->GetNumberOfComponents());
        std::cout << "copy attribute " << n << " is called " << att->GetName() << '\n';
        std::cout << "It contains " << att->GetNumberOfTuples() << " tuples and each tuple contains " << att->GetNumberOfComponents() << " coefficients\n";
        for (long i=0; i<new2old.size(); ++i) {
            att->InsertNextTuple(data->GetTuple(new2old[i]));
        }
        newpointdata->AddArray(att);
    }

    // Pass along cell attributes (unchanged)
    VTK_SMART(vtkCellData) celldata = input->GetCellData();
    std::cout << "there are " << natts << " cells attributes\n";
    natts = celldata->GetNumberOfArrays();
    for (int n=0; n<natts; ++n) {
        VTK_SMART(vtkAbstractArray) att = celldata->GetAbstractArray(n);
        std::cout << "attribute #" << n << " is called " << att->GetName() << '\n';
        std::cout << "this data set contains " << att->GetNumberOfTuples() << " tuples and each tuple contains " << att->GetNumberOfComponents() << " coefficients\n";
        shrunk->GetCellData()->AddArray(att);
    }
    return shrunk;
}

template<typename T=double>
std::ostream &operator<<(std::ostream &os, const typename std::vector<T>::iterator &i) {
    os << &i;
    return os;
}

int main(int argc, const char* argv[]) 
{
    cxxopts::Options options("filter_crease", "Manipulate crease surface");
    options.add_options()
        ("i,input", "Input filename", cxxopts::value<std::string>())
        ("o,output", "Output filename", cxxopts::value<std::string>())
        // no default values
        ("s,size", "Lower bound on CC size", cxxopts::value<int>())
        ("n,number", "Max number of CCs", cxxopts::value<int>())
        ("value", "Min scalar value", cxxopts::value<std::string>())
        ("strength", "Min ridge strength (>0)", cxxopts::value<std::string>())
        ("namev", "value name", cxxopts::value<std::string>()->default_value("values"))
        ("names", "strength name", cxxopts::value<std::string>()->default_value("ridge_strength"))
        ("pre", "Apply value filters before CC filters", cxxopts::value<bool>())
        ("t,translate", "Translate mesh", cxxopts::value<std::array<double,3>>())
        ("shrink", "Shrink dataset by removing unused vertices", cxxopts::value<bool>())
        ("stats", "Compute dataset stats and exit", cxxopts::value<bool>())
        ("v,verbose", "Verbose output", cxxopts::value<bool>())
        ("h,help", "Print usage information");
    
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("input") || (!result.count("output") && !result.count("stats"))) {
        std::cout << options.help() << '\n';
        exit(0);
    }

    std::string input = result["input"].as<std::string>();
    std::string output;
    if (result.count("output")) output = result["output"].as<std::string>();

    int minsize = 0;
    int maxsize = very_large;
    int minnb = 0;
    int maxnb = very_large;
    double minstr = minus_infinity;
    double maxstr = plus_infinity;
    double minval = minus_infinity;
    double maxval = plus_infinity;

    if (result.count("size")) minsize = result["size"].as<int>();
    if (result.count("number")) maxnb = result["number"].as<int>();;

    vec3 t(0,0,0);
    if (result.count("translate")) {
        std::array<double, 3> _t = result["translate"].as<std::array<double, 3>>();
        t[0] = _t[0];
        t[1] = _t[1];
        t[2] = _t[2];
    }

    bool verbose=false;
    if (result.count("verbose")) verbose = true;

    bool doshrink = false;
    if (result.count("shrink")) doshrink = true;

    if (verbose) std::cout << "loading dataset... " << std::flush; 

    VTK_CREATE(vtkXMLPolyDataReader, reader);
    reader->SetFileName(input.c_str());
    reader->Update();
    VTK_SMART(vtkPolyData) data = reader->GetOutput();
    if (verbose) std::cout << "done.\n";

    if (result.count("stats")) {
        int natts = data->GetPointData()->GetNumberOfArrays();
        std::cout << "There are " << natts << " data arrays in input\n";
        for (int n=0; n<natts; ++n) {
            VTK_SMART(vtkAbstractArray) arr = data->GetPointData()->GetAbstractArray(n);
            // std::vector<double> vals = array2vector(arr);
            std::vector<double> vals;
            array2vector2(vals, arr);
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _mode = spurt::mode(vals.begin(), vals.end(), 1024);
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _min = spurt::min(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _max = spurt::max(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _median = spurt::median(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            std::pair<double, double> _meanvar = spurt::meanvariance(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            std::vector<double> per = spurt::percentiles(vals.begin(), vals.end(), 21);
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            std::cout << "Stats for " << arr->GetName() << ":\n";
            std::cout << "min: " << _min << "\nmax: " << _max << "\nmean: " << _meanvar.first << "\nvariance: " << _meanvar.second << "\nmedian: " << _median << "\nmode: " << _mode << '\n'; 
            std::cout << "percentiles:\n";
            std::copy(per.begin(), per.end(), std::ostream_iterator<double>(std::cout, ", "));
            std::cout << '\n';
        }
        return 0;
    }

    if (result.count("value")) {
        std::string valasstr = result["value"].as<std::string>();
        if (valasstr.back() == '%') {
            int pct = std::stoi(valasstr.substr(0, valasstr.size()-1));
            std::vector<double> vals = array2vector(data->GetPointData()->GetAbstractArray(result["namev"].as<std::string>().c_str()));
            std::sort(vals.begin(), vals.end());
            minval = vals[std::floor(pct*vals.size()/100.)];
        }
        else minval = std::stof(result["value"].as<std::string>());
    }
    if (result.count("strength")) {
        std::string strasstr = result["strength"].as<std::string>();
        if (strasstr.back() == '%') {
            int pct = 100-std::stoi(strasstr.substr(0, strasstr.size()-1));
            std::vector<double> vals = array2vector(data->GetPointData()->GetAbstractArray(result["names"].as<std::string>().c_str()));
            std::sort(vals.begin(), vals.end());
            maxstr = vals[std::floor(pct*vals.size()/100.)];
        }
        else maxstr = std::stof(result["strength"].as<std::string>());
    }

    std::cout << "thresholds are set to: value: " << minval << ", ridge strength: " << maxstr << '\n';

    if (verbose)
        std::cout << "Initially, there are " 
                  << data->GetNumberOfCells() << " cells\n";

    if (result["pre"].as<bool>()) {
        if (verbose) std::cout << "pre is TRUE\n";
        // first filter by value ...
        data = filter_by_value<double>(data, result["namev"].as<std::string>(), minval, maxval, false);

        if (verbose)
            std::cout << "After value filtering (>" << minval << "), there are " << data->GetNumberOfCells() << " cells\n";

        // ... then filter by ridge strength ...
        data = filter_by_value<double>(data, result["names"].as<std::string>(), minstr, maxstr, false);

        if (verbose) 
            std::cout << "Aftering strength filtering (<" << minstr << "), there are " << data->GetNumberOfCells() << " cells\n";

        spurt::compute_cc_sizes(data);

        if (verbose)
            std::cout << "There are initially " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";

        // ... then filter by CC size ...
        data = filter_by_value<int>(data, "CCsizes", minsize, maxsize, true);

        if (verbose) 
            std::cout << "After filtering by CC size (>" << minsize << "), there are " << data->GetNumberOfCells() << " cells and " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";

        // ... then filter by number of CCs
        data = filter_by_value<int>(data, "CCIDs", minnb, maxnb-1, true);

        if (verbose) 
            std::cout << "After filtering by number of CCs (<" << maxnb << "), there are " << data->GetNumberOfCells() << " cells\n";
    }
    else {
        if (verbose) std::cout << "pre is FALSE\n";
        if (verbose) std::cout << "computing connected components and their sizes\n";
        spurt::compute_cc_sizes(data);

        if (verbose) {
            std::cout << "There are initially " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";
            std::cout << "filtering by size...\n";
            std::cout << "filtering size range is " << minsize << " to " << maxsize << '\n';
        }

        // First filter by CC size ...
        data = filter_by_value<int>(data, "CCsizes", minsize, maxsize, true);

        if (verbose)
            std::cout << "After filtering by CC size (>" << minsize << "), there are " << data->GetNumberOfCells() << " cells and " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";

        // ... then filter by number of CCs ...
        if (verbose) std::cout << "filtering by number of connected components, between " << minnb << " and " << maxnb << '\n';
        data = filter_by_value<int>(data, "CCIDs", minnb, maxnb-1, true);

        if (verbose) 
            std::cout << "After filtering by number of CCs (<" << maxnb << "), there are " << data->GetNumberOfCells() << " cells\n";

        // ... then filter by value ...
        data = filter_by_value<double>(data, result["namev"].as<std::string>(), minval, maxval, false);

        if (verbose)
            std::cout << "After value filtering (>" << minval << "), there are " << data->GetNumberOfCells() << " cells\n";

        // ... then filter by ridge strength.
        data = filter_by_value<double>(data, result["names"].as<std::string>(), minstr, maxstr, false);

        if (verbose)
            std::cout << "Aftering strength filtering (<" << minstr << "), there are " << data->GetNumberOfCells() << " cells\n";
    }

    if (doshrink) {
        if (verbose) std::cout << "shrinking dataset...\n";
        data = shrink(data);
    }

    if (verbose)
        std::cout << "exporting filtered crease mesh in " << output << "... " << std::flush;
    VTK_CREATE(vtkXMLPolyDataWriter, writer);
    writer->SetFileName(output.c_str());
    writer->SetInputData(data);
    writer->Write();
    if (verbose) std::cout << "done.\n";
    return 0;
}

