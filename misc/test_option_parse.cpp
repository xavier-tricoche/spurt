#include <array>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <misc/option_parse.hpp>
#include <math/fixed_vector.hpp>


struct time_of_day {
    int hour, minute, second;
};

std::ostream& operator<<(std::ostream& out, time_of_day& tod)
{
    out << tod.hour << ':' << tod.minute << ':' <<  tod.second;
    return out;
}

void validate(time_of_day& tod, const std::vector<std::string>& s) {
    // Expected format:
    // X:Y:Z
    
    namespace po = boost::program_options;
    
    std::cout << "****\n****\n****\nentering custom validate:\n";
    for (int i=0; i<s.size(); ++i) {
        std::cout << "parameter #" << i+1 << ": " << s[i] << '\n';
    }
    
    // tokenize input string
    std::vector<std::string> tokenized;
    for (int i=0; i<s.size(); ++i) {
        char * cstr = new char [s[i].length()+1];
        std::strcpy (cstr, s[i].c_str());
        char* token=std::strtok(cstr, ": ");
        while (token!=NULL) {
            std::cout << "token=" << token << '\n';
            tokenized.push_back(token);
            token=std::strtok(NULL, ": ");
        }
    }
    if (tokenized.empty()) {
        throw po::invalid_syntax(po::invalid_syntax::missing_parameter);
    }
    else if (tokenized.size()>3) {
        throw po::invalid_syntax(po::invalid_syntax::extra_parameter);
    }
    tod.hour = std::stoi(tokenized[0]);
    if (tokenized.size()>1) {
        tod.minute = std::stoi(tokenized[1]);
        if (tokenized.size()==3) {
            tod.second = std::stoi(tokenized[2]);
        }
    }
}

// typedef spurt::command_line::custom_parser< time_of_day >::parser_type parser_t;

// void validate(boost::any& v, const std::vector<std::string>& values,
//               time_of_day* target_type, int) {
//
//     std::cout << "****\n****\n****\nEntering time_of_day validate with "
//         << values.size() << " strings\n";
//
//     if (v.empty()) {
//         v = boost::any(time_of_day());
//     }
//     time_of_day* tv = boost::any_cast<time_of_day>(&v);
//     assert(NULL != tv);
//     validate(*tv, values);
// }

void validate(boost::any& v, const std::vector<std::string>& values,
              time_of_day* target_type, long) {
                  
    std::cout << "****\n****\n****\nEntering time_of_day validate with "
        << values.size() << " strings\n";
    
    if (v.empty()) {
        v = boost::any(time_of_day());
    }
    time_of_day* tv = boost::any_cast<time_of_day>(&v);
    assert(NULL != tv);
    validate(*tv, values);
}

namespace boost { namespace program_options {

// template<>
// void validate<>(boost::any& v, const std::vector<std::string>& values,
//                 time_of_day* target_type, int) {
//
//     std::cout << "****\n****\n****\nEntering time_of_day validate within boost::po namespace with "
//         << values.size() << " strings\n";
//
//     if (v.empty()) {
//         v = boost::any(time_of_day());
//     }
//     time_of_day* tv = boost::any_cast<time_of_day>(&v);
//     assert(NULL != tv);
//     validate(*tv, values);
// }

template<>
void validate<>(boost::any& v, const std::vector<std::string>& values,
                time_of_day* target_type, long) {
              
    std::cout << "****\n****\n****\nEntering time_of_day validate within boost::po namespace with "
        << values.size() << " strings\n";

    if (v.empty()) {
        v = boost::any(time_of_day());
    }
    time_of_day* tv = boost::any_cast<time_of_day>(&v);
    assert(NULL != tv);
    validate(*tv, values);
}
}

}

int main(int argc, const char* argv[]) {

    namespace xcl = spurt::command_line;
    
    char c = 'a';
    short sh = 1;
    int i = 1;
    float f = 1.1;
    double d = 0.1;
    double dd = 0.13456789;
    
    // parser_t validator = validate;
    
    std::string str2, str, str3 = "abcd";
    nvis::fixed_vector<int, 3> iv(1, 2, 3);
    nvis::fvec2 fv(1.1, 2.2);
    nvis::vec4 dv(0.1, 0.2, 0.3, 0.4);
    std::array<double, 5> da = {6, 5, 4, 3, 2};
    std::vector<short> svec;
    std::vector<std::string> tod_str;
    typedef std::vector<short>::iterator iter;
    time_of_day t_o_d;
    
    xcl::option_traits 
        required_group(true, false, "Required Options"), 
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");
        
    xcl::option_parser parser(argv[0],
        "Testing spurt::commandline::option_parser");
        
    try {
        // parser.use_default_symbol();
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("letter", c, 'b', "A character");
        parser.add_value("word", str3, str3, "A string", optional_group);
        parser.add_value("number", sh, 3, "A short");
        parser.add_value("input", str, "Input file (a string)", positional_group);
        parser.add_value("output", str2, "Output file (a string)", required_group);
        parser.add_value("other_number", i, 4, "Some random number");
        parser.add_value("real", f, 2.2, "Some floating point number");
        parser.add_value("double", dd, dd, "Some double floating point number");
        parser.add_sequence("numbers", svec, "A series of numbers");
        parser.add_tuple<3>("3tuple", iv, iv, "A tuple of 3 integers");
        parser.add_tuple<2>("2tuple", fv, fv, "A tuple of 2 floats");
        parser.add_tuple<4>("4tuple", dv, dv, "A tuple of 4 doubles");
        // parser.add_tuple<5>("5tuple", da, da, "A tuple of 5 doubles");
        // parser.add_sequence("tod", tod_str, "A time of day: X:Y:Z", optional_group, "<time>");
        // parser.add_custom< time_of_day >("tod", t_o_d, "A custom parameter type: time of day: X:Y:Z", 1, required_group, "time");
        
        parser.parse(argc, argv);
        // validate(t_o_d, tod_str);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options enteredso far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    std::cout << "Parameter values:\n"
        << "c     = " << c << '\n'
        << "sh    = " << sh << '\n'
        << "i     = " << i << '\n'
        << "f     = " << f << '\n'
        << "d     = " << d << '\n'
        << "st    = " << str << '\n'
        << "st2   = " << str2 << '\n'
        << "iv    = " << iv << '\n'
        << "fv    = " << fv << '\n'
        << "dv    = " << dv << '\n'
        << "tod   = " << t_o_d << '\n';
    std::for_each(svec.begin(), svec.end(), [] (short s) { std::cout << s << ", "; });
    std::cout << "\n";
    
    return 0;
}
