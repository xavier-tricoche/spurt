#include <iostream>
#include <iomanip>
#include <string>
#include <image/nrrd_wrapper.hpp>
#include <misc/strings.hpp>

int main(int argc, char* argv[]) {
    std::string filename = argv[1];
    
    Nrrd* nin = spurt::readNrrd(filename);
    int nbcomments = nin->cmtArr->len;
    std::cout << nbcomments << " comments in input\n";
    std::cout << std::setprecision(12);
    for (int i=0; i<nbcomments; ++i) {
        std::cout << "comment #" << i << ": " << nin->cmt[i] << '\n';
        std::vector<std::string> strs;
        spurt::tokenize(strs, nin->cmt[i], " =[],");
        if (strs[0] == "input") {
            std::cout << "input file = " << strs[2] << '\n';
        }
        else if (strs[0] == "initial") {
            std::cout << "initial time = " << std::stod(strs[2]) << '\n';
        }
        else if (strs[0] == "current") {
            std::cout << "current time = " << std::stod(strs[2]) << '\n'; 
        }
        else if (strs[0] == "epsilon") {
            std::cout << "epsilon = " << std::stod(strs[1]) << '\n';
        }
        else if (strs[0] == "resolution") {
            std::vector<std::string> _strs;
            spurt::tokenize(_strs, strs[1], "x");
            std::cout << "resolution = " << std::stoi(_strs[0]) << " by "
                << std::stoi(_strs[1]) << '\n';
        }
        else if (strs[0] == "bounds") {
            std::vector<std::string> _strs;
            std::cout << "bounds are: longitude from " << std::stod(strs[1])
                      << " to " << std::stod(strs[4]) << ", latitude from "
                      << std::stod(strs[2]) << " to " << std::stod(strs[5])
                      << "\n";
        }
        else if (strs[0] == "dt") {
            std::cout << "dt = " << std::stod(strs[1]) << '\n';
        }
    }
    nrrdNuke(nin);
    
    return 1;
}

