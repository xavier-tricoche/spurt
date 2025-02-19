#include <boost/filesystem.hpp>
#include <string>
#include <vector>
#include <regex>
#include <iostream>

namespace fs = boost::filesystem;

std::vector<std::string> 
regex_match_files(const std::string& basename) {
    fs::path apath(basename);
    fs::path directory = apath.parent_path();
    fs::path search_path = apath.filename().string();
    std::string search_str = search_path.string();

    std::vector<std::string> matching_files;
    std::regex pattern(search_str);

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (fs::is_regular_file(entry)) {
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, pattern)) {
                matching_files.push_back(entry.path().string());
            }
        }
    }
    return matching_files;
}

int main(int argc, char** argv) {
    std::vector<std::string> matching = regex_match_files(argv[1]);
    std::cout << "matching files\n";
    for (auto str : matching) { std::cout << str << '\n'; }
}