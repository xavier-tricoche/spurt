#include <iostream>
#include <fstream>
#include <math/types.hpp>
#include <vector>
#include <list>
#include <teem/nrrd.h>


std::vector<spurt::fvec3> points;
std::vector<spurt::fvec3> colors;
std::vector<std::vector<unsigned int> > lines;

char* input, *output;
float sat;
int minsize;
float col[3];
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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &input,       NULL,       "input VTK file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  0, 1, &output,      "none",     "output VTK file name");
    hestOptAdd(&hopt, "min",    "min size",         airTypeInt,     0, 1, &minsize,     "0",        "threshold on line size");
    hestOptAdd(&hopt, "s",      "saturation",       airTypeFloat,   0, 1, &sat,         "1",        "color saturation");
    hestOptAdd(&hopt, "c",      "fixed color",      airTypeFloat,   3, 3, col,          "-1 -1 -1", "fixed color to be assigned to all curves");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Assign colors to points lying on line objects based on their orientation",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

spurt::fvec3 color(unsigned int lid, unsigned int pid)
{
    const std::vector<unsigned int>& l = lines[lid];
    spurt::fvec3 dir;
    if (l.size() < 2) {
        return spurt::fvec3(0, 0, 0);
    }
    if (pid == 0) {
        dir = points[l[pid+1]] - points[l[pid]];
    } else if (pid == lines[lid].size() - 1) {
        dir = points[l[pid]] - points[l[pid-1]];
    } else {
        dir = points[l[pid+1]] - points[l[pid-1]];
    }
    
    if (spurt::norm(dir) == 0) {
        return spurt::fvec3(0, 0, 0);
    }
    return spurt::abs(dir / spurt::norm(dir));
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    std::string in(input), out(input);
    if (strcmp(output, "none")) {
        out = std::string(output);
    }
    
    std::fstream vtk_in(in.c_str(), std::ios::in);
    std::fstream vtk_out(out.c_str(), std::ios::out);
    
    // parse header
    std::string buffer, dummy;
    for (int i = 0 ; i < 4 ; ++i) {
        std::getline(vtk_in, buffer);
        vtk_out << buffer << '\n';
    }
    // parse points
    int npts;
    vtk_in >> buffer >> npts >> dummy;
    vtk_out << buffer << " " << npts << " " << dummy << '\n';
    std::getline(vtk_in, buffer);
    points.resize(npts);
    colors.resize(npts);
    std::fill(colors.begin(), colors.end(), spurt::fvec3(0, 0, 0));
    for (int i = 0 ; i < npts ; ++i) {
        spurt::fvec3& p = points[i];
        vtk_in >> p[0] >> p[1] >> p[2];
        vtk_out << p[0] << " " << p[1] << " " << p[2] << '\n';
    }
    std::getline(vtk_in, buffer);
    // parse lines
    int nlines, n;
    vtk_in >> buffer >> nlines >> dummy;
    vtk_out << buffer << " " << nlines << " " << dummy << '\n';
    std::getline(vtk_in, buffer);
    lines.resize(nlines);
    for (int i = 0 ; i < nlines ; ++i) {
        vtk_in >> n;
        std::vector<unsigned int>& l = lines[i];
        l.resize(n);
        for (int j = 0 ; j < n ; ++j) {
            vtk_in >> l[j];
        }
        
        if (n >= minsize) {
            vtk_out << n;
            for (int j = 0 ; j < n ; ++j) {
                vtk_out << " " << lines[i][j];
                colors[lines[i][j]] = color(i, j);
            }
            vtk_out << '\n';
        }
    }
    vtk_in.close();
    
    vtk_out << "POINT_DATA " << npts << '\n';
    vtk_out << "COLOR_SCALARS " << npts << " 3\n";
    
    bool have_color = (*std::min_element(&col[0], &col[3]) >= 0);
    for (int i = 0 ; i < colors.size() ; ++i) {
        if (!have_color) {
            vtk_out << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << '\n';
        } else {
            vtk_out << col[0] << " " << col[1] << " " << col[2] << '\n';
        }
    }
    vtk_out.close();
}
















