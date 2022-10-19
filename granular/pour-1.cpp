#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <iomanip>
#include <data/locator.hpp>
#include <math/fixed_vector.hpp>
#include <limits>
#include <iostream>

std::string out_name, format;
double minx=0, miny=0, minz=0, maxx=0.24, maxz=0.24, maxy=1.0, r=0.01, d=1180;
int nb=3456;
bool periodic = true;

std::string to_fortran(double x) {
	int e = 0;
	float m(x);
	if (x>1) for (; m>1.0 ; ++e) m/=10.;
	else for (; m<0.1 ; --e) m*=10;
	std::ostringstream os;
	os << m;
	while (os.str().size() < 8) {
		os << '0';
	}
	os << "E" << (e >= 0 ? '+' : '-') << std::setw(2) << std::setfill('0') << abs(e);
	return os.str();
}

void printUsageAndExit( const std::string& argv0, const std::string& offending="", 
						bool doExit = true )
{
	if (offending != "") {
		std::cerr << "ERROR: " << offending << std::endl;
	}
  	std::cerr 
    << "Usage  : " << argv0 << " [parameters] [options]\n"
	<< "Parameters:\n"
	<< "    -o  | --output <basename>         Output file base name\n"
	<< "Options:\n"
	<< "    -h  | --help                      Print this information\n"
	<< "    -r  | --radius <float>            Particle radius\n"
	<< "    -n  | --number <int>              Number of particles in system\n"
	<< "    -bx | --boundsx <float> <float>   Domain bounding box along x-axis (horizontal)\n"
	<< "    -by | --boundsy <float> <float>   Domain bounding box along y-axis (horizontal)\n"
	<< "    -bz | --boundsz <float> <float>   Domain bounding box along z-axis (vertical)\n"
	<< "    -p  | --periodic <bool>           Domain is periodic in horizontal plane\n"
	<< "    -f  | --format <string>           Export format (\"MDS\" or \"LIGGGHTS\" or \"BOTH\")\n"
	<< "    -d  | --density <float>           Particle material density\n"
	<< std::endl;
	
	if (doExit) exit(1);
}

int main(int argc, char* argv[])
{
	out_name = "none specified";
	format = "BOTH";
	for (int i=1; i<argc ; ++i) {
		std::string arg(argv[i]);
		if (arg == "-o" || arg == "--output") {
			if (i == argc-1) {
				printUsageAndExit(argv[0], "missing output");
			}
			out_name = argv[++i];
		}		
		else if (arg == "-h" || arg == "--help") {
			printUsageAndExit(argv[0]);
		}
		else if (arg == "-n" || arg == "--number") {
			if (i == argc-1) {
				printUsageAndExit(argv[0], "missing number");
			}
			nb = atoi(argv[++i]);
		}
		else if (arg == "-bx" || arg == "--boundx") {
			if (i >= argc-2) {
				printUsageAndExit(argv[0], "missing X bounds");
			}
			minx = atof(argv[++i]);
			maxx = atof(argv[++i]);
		}	
		else if (arg == "-by" || arg == "--boundy") {
			if (i >= argc-2) {
				printUsageAndExit(argv[0], "missing Y bounds");
			}
			miny = atof(argv[++i]);
			maxy = atof(argv[++i]);
		}	
		else if (arg == "-bz" || arg == "--boundz") {
			if (i >= argc-2) {
				printUsageAndExit(argv[0], "missing Z bounds");
			}
			minz = atof(argv[++i]);
			maxz = atof(argv[++i]);
		}
		else if (arg == "-p" || arg == "--periodic") {
			if (i == argc-1) {
				printUsageAndExit(argv[0], "missing period boolean");
			}
			periodic = atoi(argv[++i]);
		}		
		else if (arg == "-r" || arg == "--radius") {
			if (i == argc-1) {
				printUsageAndExit(argv[0], "missing radius");
			}
			r = atof(argv[++i]);
		}		
		else if (arg == "-f" || arg == "--format") {
			if (i == argc-1) {
				printUsageAndExit(argv[0], "missing format");
			}
			format = argv[++i];
			if (format != "NJIT" && format != "MDS" && format != "LIGGGHTS" && format != "BOTH") {
				printUsageAndExit(argv[0], "unrecognized output format");
			}
		}
		else if (arg == "-d" || arg == "--density") {
			if (i == argc-1) {
				printUsageAndExit(argv[0], "missing density");
			}
			d = atof(argv[++i]);
		}
		else {
		       printUsageAndExit(argv[0], "invalid argument");
		}
	}
	if (out_name == "none specified") {
		printUsageAndExit(argv[0], "missing output file name");
	}
	double delta_x = maxx-minx;
	double delta_y = maxy-miny;
	double delta_z = maxz-minz;
	
	typedef xavier::point_locator<double, int, 3> 	locator_type;
	typedef locator_type::point_type				point_type;
	
	locator_type locator;
	std::vector<nvis::vec3> points;
	
	srand48(time(0));
	while (points.size() < nb) {
		nvis::vec3 p(minx + drand48()*delta_x,
					miny + drand48()*delta_y,
					minz + drand48()*delta_z);
		if (points.size()) {
			point_type pt = locator.find_nearest_point(p);
			if (nvis::norm(p-pt.coordinate()) <= 2.*r) continue;
		}	
			
		points.push_back(p);
		locator.insert(point_type(p, 0));
		if (periodic) {
			int closex = 0, closey = 0;
			if (p[0] <= minx+r) {
				locator.insert(point_type(p+nvis::vec3(delta_x, 0, 0), 0));
				closex = 1;
			}
			else if (p[0] >= maxx-r) {
				locator.insert(point_type(p-nvis::vec3(delta_x, 0, 0), 0));
				closex = 2;
			}
			if (p[1] <= miny-r) {
				locator.insert(point_type(p+nvis::vec3(0, delta_y, 0), 0));
				closey = 1;
			}
			else if (p[1] >= maxy-r) {
				locator.insert(point_type(p-nvis::vec3(0, delta_y, 0), 0));
				closey = 2;
			}
			if (closex == 1 && closey == 1) {
				locator.insert(point_type(p+nvis::vec3(delta_x, delta_y, 0), 0));
			}
			else if (closex == 1 && closey == 2) {
				locator.insert(point_type(p+nvis::vec3(delta_x, -delta_y, 0), 0));
			}
			else if (closex == 2 && closey == 1) {
				locator.insert(point_type(p+nvis::vec3(-delta_x, delta_y, 0), 0));
			}
			else if (closex == 2 && closey == 2) {
				locator.insert(point_type(p+nvis::vec3(-delta_x, -delta_y, 0), 0));
			}
		}
	}
	
	std::ostringstream os;
	os << out_name << "_mds.txt";
	if (format == "NJIT" || format == "MDS" | format == "BOTH") {
		std::fstream out(os.str().c_str(), std::ios::out);
		for (int i=0 ; i<points.size() ; ++i) {
			const nvis::vec3& p = points[i];
			out << std::setw(5) << std::setfill(' ') << i << "  " 
			<< to_fortran(p[0]) << "  " 
			<< to_fortran(p[2]) << "  " 
			<< to_fortran(p[1]) << '\n';
		}
		out.close();
	}
	if (format == "LIGGGHTS" || format == "BOTH") {
		os.clear();
		os.str("");
		os << out_name << "_liggghts.txt";
		std::fstream out(os.str().c_str(), std::ios::out);
		out << "# Atom file created by " << argv[0] << "\n\n";
		out << nb << " atoms\n\n";
		out << "1 atom types\n\n";
		out << minx << " " << maxx << " xlo xhi\n";
		out << miny << " " << maxy << " ylo yhi\n";
		out << minz << " " << maxz << " zlo zhi\n\n";
		out << "Atoms\n\n";
		for (int i=0 ; i<points.size() ; ++i) {
			const nvis::vec3& p = points[i];
			out << i+1 << " 1 " << 2*r << " " << d << " " << p[0] << " " << p[1] << " " << p[2] << '\n';
		}
		out.close();
	}
}