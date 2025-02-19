#include <vtk/vtk_utils.hpp>
#include "filter_ccs.hpp"
#include <misc/cxxopts.hpp>
#include <string>

int main(int argc, const char* argv[]) 
{
    cxxopts::Options options("test_filter_ccs", "Test mesh pruning by size of connected components");
    options.add_options()
        ("i,input", "Input filename", cxxopts::value<std::string>())
        ("o,output", "Output filename", cxxopts::value<std::string>())
        ("s,size", "Lower bound on CC size", cxxopts::value<int>()->default_value("1000"))
        ("n,number", "Max number of CCs", cxxopts::value<int>()->default_value("10"))
        ("h,help", "Print usage information");
    
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("input") || !result.count("output")) {
        std::cout << options.help() << '\n';
        exit(0);
    }

    std::string input = result["input"].as<std::string>();
    std::string output = result["output"].as<std::string>();
    int minsize = result["size"].as<int>();
    int maxnb = result["number"].as<int>();

    std::cout << "loading dataset... " << std::flush; 
    VTK_CREATE(vtkXMLPolyDataReader, reader);
    reader->SetFileName(input.c_str());
    reader->Update();
    VTK_SMART(vtkPolyData) data = reader->GetOutput();
    std::cout << "done.\n";

    std::cout << "pruning mesh with min size = " << minsize << "... " << std::flush;
    spurt::prune_mesh(data, minsize, maxnb);
    std::cout << "done.\n";

    std::cout << "exporting pruned mesh... " << std::flush;
    VTK_CREATE(vtkXMLPolyDataWriter, writer);
    writer->SetFileName(output.c_str());
    writer->SetInputData(data);
    writer->Write();
    std::cout << "done.\n";
    return 0;
}

