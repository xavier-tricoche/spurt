#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include <VTK/vtk_utils.hpp>
#include <VTK/vtk_camera_helper.hpp>
#include <misc/option_parse.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <format/filename.hpp>
#include <math/stat.hpp>

#include <vtkDataSet.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>
#include <vtkPlanes.h>
#include <vtkCommand.h>
#include <vtkImageData.h>
#include <vtkNrrdReader.h>
#include <vtkContourFilter.h>

typedef nvis::fixed_vector<double, 3> color_type;
typedef nvis::fixed_vector<double, 3> alpha_region;

std::string input_name, alpha_name, color_name, img_name, cam_in, cam_out, csv_name;
bool verbose=false;
color_type bg_col(0, 0, 0);
nvis::ivec2 res(1280, 800);
nvis::ivec3 dims;
double dist=0.1;
double stddev=0;
double signal_min=5;
double feature_min=7;
double volume;
std::vector<double> alpha_tf;
std::vector<double> color_tf;
vtkDoubleArray* scalars;
std::string input_alpha, input_color;
nvis::bbox3 global_bounds, whole_bounds, crown_bounds, germ_bounds, wall_bounds;
vtkSmartPointer<vtkImageData> image, filtered;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> window;
vtkSmartPointer<vtkColorTransferFunction> colors;
nvis::ivec3 minc, maxc;
std::vector< std::pair<double, double> > alpha_cps;
std::vector< std::pair<double, color_type > > color_cps;
std::string whole_string, germ_string, crown_string, wall_string;

void read_csv(const std::string& csv_filename, int soaktime, int grain) {
	std::fstream input(csv_filename, std::ios::in);
	if (!input) {
		std::cout << "unable to open " << csv_filename << '\n';
		exit(1);
	}
	int col = soaktime/6*3 + grain;
	std::string aline;
	for (int i=0; i<4; ++i) {
		std::getline(input, aline); // <what>,sum ...
		std::cout << "Read: " << aline << '\n';
		size_t nchars = aline.find_first_of(",");
		std::string what = aline.substr(0, nchars);
		std::cout << "what=" << what << '\n';
		std::getline(input, aline); // ,mean...
		std::getline(input, aline); // ,variance...
		std::getline(input, aline); // ,bounds
		std::cout << "bounds line: " << aline << '\n';
		size_t s=0;
		for (int c=0; c<col; ++c) { // skip as many bounds strings as necessary
			s=aline.find_first_of("\"", ++s); std::cout << "s=" << s << '\n';
			s=aline.find_first_of("\"", ++s); std::cout << "s=" << s << '\n';
		}
		size_t begin = aline.find_first_of("\"", ++s);
		size_t end = aline.find_first_of("\"", begin+1);
		std::cout << "begin=" << begin << ",end=" << end << '\n';
 		std::string bounds_str = aline.substr(begin+1, end-begin-1);;
		if (what == "whole" ) {
			whole_string = bounds_str;
			if (verbose) std::cout << "whole=\"" << whole_string << "\"\n";
		}
		else if (what == "Crown") {
			crown_string = bounds_str;
			if (verbose) std::cout << "crown=\"" << crown_string << "\"\n";
		}
		else if (what == "Germ") {
			germ_string = bounds_str;
			if (verbose) std::cout << "germ=\"" << germ_string << "\"\n";
		}
		else if (what == "Vitreous + germ" || what == "Vitreous + Germ") {
			wall_string = bounds_str;
			if (verbose) std::cout << "wall=\"" << wall_string << "\"\n";
		}
		else {
			std::cout << "unrecognized item: " << what << '\n';
			exit(1);
		}
	}
}


inline nvis::vec3 vec3(double* v) {
	return nvis::vec3(v[0], v[1], v[2]);
}

inline nvis::ivec3 ivec3(int* v) {
    return nvis::ivec3(v[0], v[1], v[2]);
}

inline nvis::ivec3 id2c(int id, const nvis::ivec3& dims) {
    nvis::ivec3 coord;
    coord[0] = id % dims[0];
    id /= dims[0];
    coord[1] = id % dims[1];
    coord[2] = id / dims[1];
    return coord;
}

inline std::string bounds_to_str(const nvis::bbox3& b) {
	std::ostringstream os;
	os << b.min() << " -> " << b.max();
	return os.str();
}

std::ostream& bold_on(std::ostream& os) {
	return os << "\e[1m";
}
std::ostream& bold_off(std::ostream& os) {
	return os << "\e[0m";
}

nvis::bbox3 str_to_bbox(const std::string& str) {
    std::string copy(str);
    for (int i=0; i<str.size(); ++i) {
        if (str[i] == '[' || str[i] == ']' || str[i] == ',') {
            copy[i] = ' ';
        }
        else if (i<str.size()-1 && str[i] == '-' && str[i+1] == '>') {
            copy[i++] = ' ';
            copy[i] = ' ';
        }
    }
    nvis::bbox3 box;
    std::istringstream iss(copy);
    for (int i=0; i<6; ++i) {
        if (i<3) iss >> box.min()[i];
        else iss >> box.max()[i-3];
        if (iss.fail()) {
            std::cerr << "unable to parse bounding box information.\n";
            exit(1);
        }
    }
    return box;
}

void initialize(int argc, char* argv[]) {
    namespace xcl = xavier::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Partial measure and visualization of scalar volume");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", input_name, "Input filename",
		required_group);
		parser.add_value("csv", csv_name, "CSV filename", optional_group);
        parser.add_value("whole", whole_string, "whole bounds", optional_group);
        parser.add_value("crown", crown_string, "crown bounds", optional_group);
        parser.add_value("germ", germ_string, "germ bounds", optional_group);
        parser.add_value("wall", wall_string, "germ + wall region", optional_group);
        parser.add_value("alpha", alpha_name,
		"Transfer function filename or inline definition", optional_group);
        parser.add_value("color", color_name,
						 "Transfer function filename or inline definition",
						  optional_group);
        parser.add_value("min", signal_min, signal_min, "Signal/noise threshold", optional_group);
        parser.add_value("fmin", feature_min, feature_min, "Feature threshold", optional_group);
        parser.add_value("save", img_name, "Snapshot filename",
		optional_group);
        parser.add_value("camin", cam_in, "Camera input filename",
		optional_group);
        parser.add_value("camout", cam_out, "Camera output filename",
		optional_group);
        parser.add_value("bg", bg_col, bg_col, "Background color",
		optional_group);
        parser.add_tuple<2>("res", res, res, "Image resolution",
		optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output",
		optional_group);

        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception: "
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}


bool selected_region(const nvis::bbox3& bounds, nvis::ivec3& min, nvis::ivec3& max) {
	// make sampling region fit within dataset bounds
	nvis::vec3 _min = nvis::max(bounds.min(), global_bounds.min());
	nvis::vec3 _max = nvis::min(bounds.max(), global_bounds.max());
	if (nvis::any(_min > _max)) {
		std::cerr << "invalid sampling region: " << bounds_to_str(bounds) << '\n';
        std::cout << "min=" << _min << ", max=" << _max << '\n';
        std::cout << "global bounds: " << bounds_to_str(global_bounds) << '\n';
		return false;
	}
    int minid = image->FindPoint(_min[0], _min[1], _min[2]);
    int maxid = image->FindPoint(_max[0], _max[1], _max[2]);
	int tmp = minid;
	minid = std::min(minid, maxid);
	maxid = std::max(tmp, maxid);
	// std::cout << "minid=" << minid << ", maxid=" << maxid << '\n';
    min = id2c(minid, dims);
    max = id2c(maxid, dims);
	return true;
}

double stats(const std::vector<double>& values, double& min, double& max,
             double& mean, double& stddev) {
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();
	double sum=0;
	std::for_each(values.begin(), values.end(), [&](double v)
		{
			sum += v;
            if (v<min) min=v;
            if (v>max) max=v;
		});
	mean = sum/static_cast<double>(values.size());
	double variance=0;
	std::for_each(values.begin(), values.end(), [&](double v)
		{
			variance += (v-mean)*(v-mean);
		});
	variance /= static_cast<double>(values.size()-1);
	stddev = sqrt(variance);
	return sum;
}

nvis::fixed_vector<double, 5>  measure(const nvis::bbox3& bounds,
                                       double threshold) {
    if (!selected_region(bounds, minc, maxc)) {
    	return 0;
    }
	if (verbose) {
		std::cout << "overall domain: " << bounds_to_str(global_bounds) << '\n';
		std::cout << "sampled region: " << bounds_to_str(bounds)  << '\n';
		std::cout << "image region: " << minc << " -> " << maxc << '\n';
		std::cout << "signal threshold: " << signal_min << '\n';
		std::cout << "feature threshold: " << threshold << '\n';
		// std::cout << "voxel volume:   " << volume << '\n';
	}
	std::vector<double> noise_values, feature_values, non_feature_values;
    for (int k=minc[2]; k<=maxc[2]; ++k) {
        for (int j=minc[1]; j<=maxc[1]; ++j) {
            for (int i=minc[0]; i<=maxc[0]; ++i) {
				vtkIdType id = i+dims[0]*(j+dims[1]*k);
                double value = scalars->GetTuple(id)[0];
				if (value < signal_min) {
					noise_values.push_back(value);
				}
				else if (value < threshold) {
					non_feature_values.push_back(value);
				}
				else {
					feature_values.push_back(value);
				}
            }
        }
    }
	size_t nsamples = noise_values.size() + feature_values.size() + non_feature_values.size();
	if (verbose) std::cout << "sampling (" << nsamples << ") done\n";
	double min[3], max[3], mean[3], stddev[3], sum[3];
	sum[0] = stats(noise_values, min[0], max[0], mean[0], stddev[0]);
	sum[1] = stats(non_feature_values, min[1], max[1], mean[1], stddev[1]);
	sum[2] = stats(feature_values, min[2], max[2], mean[2], stddev[2]);
	if (verbose) {
		std::cout << nsamples << " total samples\n";
		std::cout << static_cast<double>(nsamples-noise_values.size())/nsamples*100. << "% signal values\n";
		std::cout << noise_values.size() << " noise samples (" << static_cast<double>(noise_values.size())/nsamples*100. << "%)\n";
		std::cout << "sum=" << sum[0] << ", min=" << min[0] << ", max=" << max[0] << ", mean=" << mean[0] << ", stddev=" << stddev[0] << '\n';
		std::cout << non_feature_values.size() << " non feature samples (" << static_cast<double>(non_feature_values.size())/nsamples*100. << "%)\n";
		std::cout << "sum=" << sum[1] << ", min=" << min[1] << ", max=" << max[1] << ", mean=" << mean[1] << ", stddev=" << stddev[1] << '\n';
		std::cout << feature_values.size() << " feature samples (" << static_cast<double>(feature_values.size())/nsamples*100. << "%)\n";
		std::cout << "sum=" << sum[2] << ", min=" << min[2] << ", max=" << max[2] << ", mean=" << mean[2] << ", stddev=" << stddev[2] << '\n';
    }
	std::cout << "following results for selected value range x region selection:\n";
	std::cout << "sum / mean / variance / bounds:\n";
	std::cout << sum[2] << '\n';
	std::cout << mean[2] << '\n';
	std::cout << stddev[2] << '\n';
	std::cout << bounds_to_str(bounds) << '\n';
		// std::cout << ntotal << " values sampled of which " << nbad << " where outside implicit region\n";
    nvis::fixed_vector<double, 5> r;
    r[0] = sum[2];
    r[1] = min[2];
    r[2] = max[2];
    r[3] = mean[2];
    r[4] = stddev[2];
    return r;
}

double pct(double a, double b) {
    return 100.*a/b;
}

void set_parameters() {
	scalars =
		vtkDoubleArray::SafeDownCast(image->GetPointData()->GetScalars());

	image->ComputeBounds();
	double* domain = image->GetBounds();
	nvis::vec3 dmin = nvis::vec3(domain[0], domain[2], domain[4]);
	nvis::vec3 dmax = nvis::vec3(domain[1], domain[3], domain[5]);
	global_bounds.min() = dmin;
	global_bounds.max() = dmax;
	if (verbose) {
		std::cout << "domain: " << bounds_to_str(global_bounds) << '\n';
	}
    dims = ivec3(image->GetDimensions());
	if (verbose) {
		std::cout << "image dimensions = " << dims << '\n';
	}
	double* spacing = image->GetSpacing();
	volume = spacing[0] * spacing[1] * spacing[2];
	// center = 0.5*(dmin + dmax);
	// dfront = dside = dtop = dbottom = 100.;
}

int main(int argc, char* argv[]) {
    initialize(argc, argv);

	VTK_CREATE(vtkNrrdReader, reader);
	reader->SetFileName(input_name.c_str());
	reader->Update();
	image = vtkImageData::SafeDownCast(reader->GetOutput());

    set_parameters();
	
	nvis::bbox3 tmp_whole_bounds, tmp_crown_bounds, tmp_germ_bounds, tmp_wall_bounds;
	if (!csv_name.empty()) {
		std::istringstream iss(csv_name);
		std::string filename;
		int soak_time, grain;
		iss >> filename >> soak_time >> grain;
		read_csv(filename, soak_time, grain);
		tmp_whole_bounds = str_to_bbox(whole_string);
		tmp_crown_bounds = str_to_bbox(crown_string);
		tmp_germ_bounds = str_to_bbox(germ_string);
		tmp_wall_bounds = str_to_bbox(wall_string);
	}
	else {
		if (whole_string.empty() || crown_string.empty() || germ_string.empty() || wall_string.empty()) {
			std::cerr << "insufficient boundary information\n";
			exit(1);
		}
    	tmp_whole_bounds = str_to_bbox(whole_string);
    	tmp_crown_bounds = str_to_bbox(crown_string);
    	tmp_germ_bounds = str_to_bbox(germ_string);
    	tmp_wall_bounds = str_to_bbox(wall_string);
	}
    if (verbose) {
        std::cout << "kernel region: " << bounds_to_str(tmp_whole_bounds) << '\n';
        std::cout << "crown region: " << bounds_to_str(tmp_crown_bounds) << '\n';
        std::cout << "germ region: " << bounds_to_str(tmp_germ_bounds) << '\n';
        std::cout << "wall region: " << bounds_to_str(tmp_wall_bounds) << '\n';
    }

    whole_bounds.min()[0] = crown_bounds.min()[0] =
        std::max(tmp_whole_bounds.min()[0], tmp_crown_bounds.min()[0]);
    crown_bounds.max()[0] = wall_bounds.min()[0] =
        std::min(tmp_crown_bounds.max()[0], tmp_wall_bounds.min()[0]);
    whole_bounds.max()[0] = wall_bounds.max()[0] =
        std::max(std::max(tmp_whole_bounds.max()[0], tmp_germ_bounds.max()[0]),
                 tmp_wall_bounds.max()[0]);
    germ_bounds.min()[0] = tmp_germ_bounds.min()[0];
    germ_bounds.max()[0] = std::min(tmp_germ_bounds.max()[0],
                                    whole_bounds.max()[0]);

    whole_bounds.min()[1] = wall_bounds.min()[1] =
        std::max(tmp_whole_bounds.min()[1], tmp_wall_bounds.min()[1]);
    whole_bounds.max()[1] = wall_bounds.max()[1] =
        std::min(tmp_whole_bounds.max()[1], tmp_wall_bounds.max()[1]);
    germ_bounds.min()[1] = tmp_germ_bounds.min()[1];
    germ_bounds.max()[1] = tmp_germ_bounds.max()[1];
    crown_bounds.min()[1] = std::max(whole_bounds.min()[1],
                                     tmp_crown_bounds.min()[1]);
    crown_bounds.max()[1] = std::min(whole_bounds.max()[1],
                                     tmp_crown_bounds.max()[1]);

    whole_bounds.min()[2] = wall_bounds.min()[2] =
        std::min(tmp_whole_bounds.min()[2], tmp_wall_bounds.min()[2]);
    whole_bounds.max()[2] = wall_bounds.max()[2] = germ_bounds.max()[2] =
        std::max(std::max(tmp_whole_bounds.max()[2], tmp_wall_bounds.max()[2]),
                 tmp_germ_bounds.max()[2]);
    crown_bounds.min()[2] = tmp_crown_bounds.min()[2];
    crown_bounds.max()[2] = tmp_crown_bounds.max()[2];
    germ_bounds.min()[2] = tmp_germ_bounds.min()[2];
    germ_bounds.max()[2] = tmp_germ_bounds.max()[2];

    std::cout << bold_on << "regions after correction:\n" << bold_off;
    std::cout << "kernel region: " << bounds_to_str(whole_bounds) << '\n';
    std::cout << "crown region: " << bounds_to_str(crown_bounds) << '\n';
    std::cout << "germ region: " << bounds_to_str(germ_bounds) << '\n';
    std::cout << "wall region: " << bounds_to_str(wall_bounds) << '\n';

    nvis::fixed_vector<double, 5> whole_vals, wall_vals, crown_vals, germ_vals;
    std::cout << bold_on << "whole\n" << bold_off;
    whole_vals = measure(whole_bounds, signal_min);

    std::cout << bold_on << "crown\n" << bold_off;
    crown_vals = measure(crown_bounds, signal_min);
    std::cout << bold_on << "germ\n" << bold_off;
    germ_vals = measure(germ_bounds, feature_min);
    std::cout << bold_on << "wall and germ\n" << bold_off;
    wall_vals = measure(wall_bounds, feature_min);

    std::cout << "Results summary: absolute values: whole / crown / wall / germ / floury\n";
    std::cout << whole_vals[0]*volume << '\n';
    std::cout << crown_vals[0]*volume << '\n';
    std::cout << (wall_vals[0] - germ_vals[0])*volume << '\n';
    std::cout << germ_vals[0]*volume << '\n';
    std::cout << (whole_vals[0] - crown_vals[0] - wall_vals[0])*volume << '\n';

    std::cout << "Results summary: relative values: crown / wall / germ / floury\n";
    std::cout << pct(crown_vals[0], whole_vals[0]) << '\n';
    std::cout << pct(wall_vals[0] - germ_vals[0], whole_vals[0]) << '\n';
    std::cout << pct(germ_vals[0], whole_vals[0]) << '\n';
    std::cout << pct(whole_vals[0] - crown_vals[0] - wall_vals[0], whole_vals[0]) << '\n';


	return 0;
}
