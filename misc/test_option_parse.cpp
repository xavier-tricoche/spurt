#include <array>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <misc/option_parse.hpp>
#include <math/types.hpp>

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
    spurt::small_vector<int, 3> iv(1, 2, 3);
    spurt::fvec2 fv(1.1, 2.2);
    spurt::vec4 dv(0.1, 0.2, 0.3, 0.4);
    std::array<double, 5> da = {6, 5, 4, 3, 2};
    std::vector<short> svec;
    std::vector<std::string> tod_str;
    typedef std::vector<short>::iterator iter;
    
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
        parser.add_tuple<5>("5tuple", da, da, "A tuple of 5 doubles");
        
        parser.parse(argc, argv);
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
        << "letter       = " << c << '\n'
        << "number       = " << sh << '\n'
        << "other_number = " << i << '\n'
        << "real         = " << f << '\n'
        << "double       = " << dd << '\n'
        << "input        = " << str << '\n'
        << "ouput        = " << str2 << '\n'
        << "3tuple       = " << iv << '\n'
        << "2tuple       = " << fv << '\n'
        << "4tuple       = " << dv << '\n'
        << "numbers      = ";
    std::for_each(svec.begin(), svec.end(), [] (short s) { std::cout << s << ", "; });
    std::cout << "\n";
    
    return 0;
}
