#ifndef __XAVIER_OPTION_PARSE_HPP__

// STL
#include <exception>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

// Boost
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

// Add validator-level support for fixed-size arrays to boost::program_options
#include <misc/option_parse_array.hpp>

namespace xavier { namespace command_line {

namespace po = boost::program_options;
    
// Simple class to store general information about a group of options
class option_traits {
public:
    option_traits(bool required=false, 
                  bool positional=false) 
        : _M_default_requirement(required), 
          _M_default_positional(positional),
          _M_group_name() {}
          
    option_traits(const std::string& group_name,
                  bool required=false, 
                  bool positional=false) 
        : _M_default_requirement(required), 
          _M_default_positional(positional),
          _M_group_name(group_name) {}
          
    void required(bool required) { _M_default_requirement = required; }
    void positional(bool positional) { _M_default_positional = positional; }
    void group_name(const std::string& name) { _M_group_name = name; }
    
    bool required() const { return _M_default_requirement; }
    bool positional() const { return _M_default_positional; }
    bool has_group() const { return !_M_group_name.empty(); }
    const std::string& group_name() const { return _M_group_name; }
    
protected:
    bool         _M_default_requirement;
    bool         _M_default_positional;
    std::string  _M_group_name;
};

// Helper class to turn a user-defined parsing function
// into a callback function compatible with boost::program_options
// post-parsing notification mechanism
template<typename _Type>
struct custom_parser {
    // callback function type expected by boost::program_options
    typedef typename boost::function1<void, std::vector<std::string>&> callback_type;
    // user-defined function performing simultaneously validation of a
    // series of strings and value assignment of provided variable
    typedef typename std::function<void (_Type&, const std::vector<std::string>&)> parser_type;
    
    static callback_type callback(_Type& variable, 
                                  parser_type parser) {
        return callback_type([&](std::vector<std::string>& s) 
                             { parser(variable, s); });
    }
};

/** option_parser: General-purpose parser for command line 
    options of arbitrary types and syntax. 
    The implementation uses boost::program_options and provides
    a number of useful extensions:
    - simplified interface
    - built-in support for fixed-size array types
    - support for arbitrary types for which the user provides a 
      strings to value conversion function.
 */
class option_parser {
    typedef po::options_description          option_group;
    typedef boost::shared_ptr<option_group>  group_ptr;
    
public:
    const std::string default_group_title;
    
    option_parser(const std::string& title,
                  const std::string& default_title="Default Options") 
        : _M_all_options(title), _M_groups(), _M_positionals(), 
          _M_current_position(0), default_group_title(default_title),
          _M_default_symbol("arg"),
          _M_use_default_symbol(false), _M_use_brackets(true), 
          _M_use_short_symbols(true) {
              
        // create a group for general options
        option_group* general = 
            new option_group(default_group_title);
        _M_groups.push_back(boost::shared_ptr<option_group>(general));
        _M_groups.back()->add_options()("help,h", "Print this message");
        _M_group_map[default_group_title] = 0;
        _M_flags_string = "-h";
    }
    
    // Parameters controlling the display of option's arguments
    bool use_default_symbol() const {
        return _M_use_default_symbol;
    }
    void use_default_symbol(bool use_it=true) {
        _M_use_default_symbol = use_it;
    }
    bool use_short_symbols() const {
        return _M_use_short_symbols;
    }
    void use_short_symbols(bool use_it=true) {
        _M_use_short_symbols = use_it;
        _M_use_default_symbol = false;
    }
    bool use_brackets() const {
        return _M_use_brackets;
    }
    void use_brackets(bool use_it=true) {
        _M_use_brackets = use_it;
    }
    const std::string& default_symbol() const {
        return _M_default_symbol;
    }
    void default_symbol(const std::string& symbol) {
        _M_default_symbol = symbol;
    }
    
    std::string usage(const std::string& filename) const;
    
    // Add a parameter
    // By default, the parameter is required and non-positional
    template<typename _Type>
    void add_value(const std::string& name,
                   _Type& variable,
                   const std::string& description,
                   const option_traits& traits=option_traits(true),
                   const std::string& symbol="");
                   
    // Add and initialize parameter
    // By default, the parameter is optional and non-positional
    template<typename _Type, typename _CompatibleType = _Type>
    void add_value(const std::string& name,
                   _Type& variable,
                   const _CompatibleType& value,
                   const std::string& description,
                   const option_traits& traits=option_traits(false),
                   const std::string& symbol="");
    
    // Add and initialize boolean flag        
    // A flag is ALWAYS optional and non-positional
    void add_flag(const std::string& name, bool& variable,
                  const std::string& description,
                  const option_traits& traits);

    // Add fixed-size random access container 
    // (e.g., boost::array / std::array)
    // By default: the tuple is required and non-positional
    template<size_t N, typename _Tuple,
             typename _Type = typename _Tuple::value_type>
    void add_tuple(const std::string& name, 
                   _Tuple& variable,
                   const std::string& description,
                   const option_traits& traits=option_traits(true),
                   const std::string& symbol="");
    
    // Add and initialize fixed-size random access container 
    // (e.g., boost::array / std::array)
    // By default: the tuple is optional and non-positional
    template<size_t N, typename _Tuple,
             typename _Type = typename _Tuple::value_type,
             typename _CompatibleTuple = _Tuple>
    void add_tuple(const std::string& name, 
                   _Tuple& variable,
                   const _CompatibleTuple& value, 
                   const std::string& description,
                   const option_traits& traits=option_traits(false),
                   const std::string& symbol="");
    
    // Add std::vector
    // By default: the vector is optional and non-positional
    template<typename _Type>
    void add_sequence(const std::string& name, 
                      std::vector<_Type>& variable, 
                      const std::string& description,
                      const option_traits& traits=option_traits(false),
                      const std::string& symbol="");
                    
    // Add parameter taking N arguments with custom validation procedure
    // By default: the parameter is optional. 
    template<typename _Type, int N>
    void add_custom(const std::string& name, 
                    _Type& variable,
                    const std::string& description,
                    typename custom_parser<_Type>::parser_type parser,
                    const option_traits& traits=option_traits(false),
                    const std::string& symbol="");
    
    void parse(int argc, char* argv[]);

private:
    typedef po::positional_options_description positionals_type;
    
    option_group                  _M_all_options;
    std::vector<group_ptr>        _M_groups;
    std::map<std::string, size_t> _M_group_map;
    positionals_type              _M_positionals;
    std::string                   _M_default_symbol;
    std::string                   _M_flags_string;
    std::string                   _M_options_string;
    std::string                   _M_positionals_string;
    int                           _M_current_position;
    bool                          _M_use_default_symbol;
    bool                          _M_use_brackets;
    bool                          _M_use_short_symbols;
    
    // returns shortest valid name of an option (with optional dash)
    std::string short_form(const std::string& name, bool with_dash=false);
    
    // returns symbol representation of option parameter(s)
    std::string parameter_string(const std::string& s, int N);
    template<typename _Type, int N>
    std::string parameter_string(const std::string& s);
    
    // add boolean switch to usage string
    void add_flag_to_usage_string(const std::string& name);
    
    // add option description to usage string
    void add_option_to_usage_string(const std::string& name, 
                                    const std::string& symbol,
                                    size_t N=1, bool required=false, 
                                    bool positional=false);
    
    // return option group associated with traits
    group_ptr get_group(const option_traits& traits);
    
    // register positional option
    void register_positional(const std::string& name, int N=1);
};

} // command_line
} // xavier

inline void xavier::command_line::option_parser::
register_positional(const std::string& name, int N) {
    if (_M_current_position == -1) {
        throw std::runtime_error("Invalid sequence of positional arguments");
    }
    if (N == -1) _M_current_position = -1;
    else _M_current_position += N;
    _M_positionals.add(name.c_str(), _M_current_position);
}

// Display grammar for option <name> with value symbol <sym>
// where <short> =
//    * "-<name>" if <name> is the one-character short name of the option
//    * "--<name>" if no short name exists for this option
// DISPLAY(<name>, <sym>) = 
//    * <FORMAT(<name>, <sym>)> + " " if required
//    * [<FORMAT(<name>, <sym>)>] + " " otherwise
// FORMAT(<name>, <sym>) = 
//    * POS(<name>, <sym>) + " | " + OPT(<name>, <sym>) if positional
//    * OPT(<name>, <sym>) otherwise
// POS(<name>, <sym>) = VAL(<sym>)
// OPT(<name>, <sym>) = "-o " + VAL(<sym>)
// VAL(<sym>) = 
//    * "<sym>" if single value
//    * "<sym> <sym> "..." <sym>" (n times) if n-tuple
//    * "<sym> ..." if open-ended sequence

inline std::string xavier::command_line::option_parser::
usage(const std::string& filename) const { 
    std::string out;
    out = filename + " " + _M_flags_string + " ";
    if (!_M_options_string.empty()) out += _M_options_string;
    if (!_M_positionals_string.empty()) out += _M_positionals_string; 
    return out;
}

// <short>
inline std::string xavier::command_line::option_parser::
short_form(const std::string& name, bool with_dash) {
    size_t n = name.find(',');
    std::string s;
    if (n != std::string::npos) {
        if (with_dash) {
            return "-" + name.substr(n-1, 1);
        }
        else {
            return name.substr(n-1, 1);
        }
    }
    else if (with_dash) {
        return "--" + name;
    }
    else {
        return name;
    }
}

// VAL(<sym>)
inline std::string xavier::command_line::option_parser::
parameter_string(const std::string& s, int N) {
    if (N < 0) {
        return s + "...";
    }
    else if (N == 1) return s;
    else if (N*s.size() + N - 1 <= 10) { // somewhat random
        std::string r = s;
        for (size_t i=1 ; i<N; ++i) {
            r += " " + s;
        }
        return r;
    }
    else {
        return s + " (x" + boost::lexical_cast<std::string>(N) + ")";
    }
}

template<typename _Type, int N>
inline std::string xavier::command_line::option_parser::
parameter_string(const std::string& s) {
    if (!s.empty()) return s;
    else if (_M_use_default_symbol) 
        return parameter_string(_M_default_symbol, N);
    else {
        return parameter_string(type_as_string<_Type>
                                    (_M_use_brackets, _M_use_short_symbols), 
                                N);
    }
}

inline void xavier::command_line::option_parser::
add_flag_to_usage_string(const std::string& name) {
    std::string s = short_form(name, false);
    if (s.size() == 1) _M_flags_string += s;
    else _M_options_string += "[" + short_form(name, true) + "] ";
}

inline void xavier::command_line::option_parser::
add_option_to_usage_string(const std::string& name, 
                           const std::string& val_str,
                           size_t N, bool required, bool positional) {
    std::string opt_str = short_form(name, true) + " " + val_str;
    std::string form_str;
    if (positional) form_str = val_str + "|" + opt_str;
    else form_str = opt_str;
    std::string usage_str;
    if (required) usage_str = form_str;
    else usage_str = "[" + form_str + "]";
    if (positional) _M_positionals_string += usage_str + " ";
    else _M_options_string += usage_str + " ";
}

inline xavier::command_line::option_parser::group_ptr 
xavier::command_line::option_parser::
get_group(const option_traits& traits) {
    const std::string& name = traits.group_name();
    if (!name.empty()) {
        std::map<std::string, size_t>::const_iterator it;
        it = _M_group_map.find(name);
        if (it != _M_group_map.end()) {
            // existing group
            return _M_groups[it->second];
        }
        else {
            // new group
            option_group* __group = 
                new option_group(name);
            _M_groups.push_back(group_ptr(__group));
            _M_group_map[name] = _M_groups.size()-1;
            return _M_groups.back();
        }
    }
    else { // default group
        return _M_groups[0];
    }
}

template<typename _Type>
inline void xavier::command_line::option_parser::
add_value(const std::string& name, _Type& variable, 
          const std::string& description,
          const option_traits& traits,
          const std::string& symbol) {
    
    const std::string s = parameter_string<_Type>(symbol);
    get_group(traits)->add_options()(
        name.c_str(), 
        get_semantic(variable, traits.required(), s),
        description.c_str()
    );
    add_option_to_usage_string(name, s, 1, traits.required(), 
                               traits.positional());
                           
   if (traits.positional()) register_positional(name);
}

template<typename _Type, typename _CompatibleType>
inline void xavier::command_line::option_parser::
add_value(const std::string& name, _Type& variable, 
          const _CompatibleType& value, const std::string& description,
          const option_traits& traits,
          const std::string& symbol) {
              
    _Type init_value = boost::lexical_cast<_Type>(value);
    const std::string s = parameter_string<_Type>(symbol);
    get_group(traits)->add_options()(
        name.c_str(), 
        get_semantic(variable, traits.required(), init_value, s),
        description.c_str()
    );
    add_option_to_usage_string(name, s, 1, traits.required(), 
                               traits.positional());
    if (traits.positional()) register_positional(name);
}

inline void xavier::command_line::option_parser::
add_flag(const std::string& name, bool& variable, 
         const std::string& description,
         const option_traits& traits) {

    variable = false;
    get_group(traits)->add_options()(
        name.c_str(), po::bool_switch(&variable), description.c_str()
    );
    add_flag_to_usage_string(name);
}

template<size_t N, typename _Tuple, typename _Type, typename _CompatibleTuple>
inline void xavier::command_line::option_parser::
add_tuple(const std::string& name, _Tuple& variable, 
          const _CompatibleTuple& value, const std::string& description,
          const option_traits& traits,
          const std::string& symbol) {
              
    _Tuple initial_value;
    for (size_t i=0 ; i<N ; ++i) {
        initial_value[i] = boost::lexical_cast<_Type>(value[i]);
    }

    std::string value_as_string = "(";
    for (size_t i=0 ; i<N-1 ; ++i) {
        value_as_string += 
            boost::lexical_cast<std::string>(value[i]) + ", ";
    }
    value_as_string += 
        boost::lexical_cast<std::string>(value[N-1]) + ")";
        
    const std::string s = parameter_string<_Type, N>(symbol);
                                  
    get_group(traits)->add_options()(
        name.c_str(), 
        get_semantic(variable, traits.required(), initial_value,
                     value_as_string, s)->multitoken(),
        description.c_str()
    );
    add_option_to_usage_string(name, s, N, traits.required(), 
                               traits.positional());
    if (traits.positional()) register_positional(name, N);
}

template<size_t N, typename _Tuple, typename _Type>
inline void xavier::command_line::option_parser::
add_tuple(const std::string& name, _Tuple& variable, 
          const std::string& description,
          const option_traits& traits,
          const std::string& symbol) {
    
    const std::string s = parameter_string<_Type, N>(symbol);
    get_group(traits)->add_options()(
        name.c_str(), 
        get_semantic(variable, traits.required(), s)->multitoken(),
        description.c_str()
    );
    add_option_to_usage_string(name, s, N, traits.required(), 
                               traits.positional());
    if (traits.positional()) register_positional(name, N);
}

template<typename _Type>
inline void xavier::command_line::option_parser::
add_sequence(const std::string& name, std::vector<_Type>& variable, 
             const std::string& description,
             const option_traits& traits,
             const std::string& symbol) {
    
    const std::string s = parameter_string<_Type, -1>(symbol);
    get_group(traits)->add_options()(
        name.c_str(), 
        get_semantic(variable, traits.required(), s)->multitoken(),
        description.c_str()
    );
    add_option_to_usage_string(name, s, -1, traits.required(), 
                               traits.positional());
    if (traits.positional()) register_positional(name, -1);
}

template<typename _Type, int N>
inline void xavier::command_line::option_parser::
add_custom(const std::string& name, _Type& variable,
           const std::string& description,
           typename custom_parser<_Type>::parser_type parser,
           const option_traits& traits,
           const std::string& symbol) {

    typedef custom_parser<_Type> parser_type;
    
    const std::string s = parameter_string<_Type>(symbol);
    get_group(traits)->add_options()(
        name.c_str(), 
        get_semantic(variable, traits.required(), s)->multitoken()
                     ->notifier(parser_type::callback(variable, parser)),
        description.c_str()
    );
    add_option_to_usage_string(name, s, N, traits.required(), 
                               traits.positional());
                               
    if (traits.positional()) register_positional(name, N);
}

inline void xavier::command_line::
option_parser::parse(int argc, char* argv[]) {
    // merge the various option groups together
    for (size_t i=0 ; i<_M_groups.size() ; ++i) {
        _M_all_options.add(*_M_groups[i]);
    }
        
    boost::filesystem::path p(argv[0]);
    std::string filename = p.filename().string();
    
    po::variables_map vm;
    try {
        po::store(
            po::command_line_parser(argc, argv).
            options(_M_all_options).
            positional(_M_positionals).run(), 
            vm
        );
                  
        if (vm.count("help")) {
            std::cout << usage(filename) << "\n\n" 
                      << _M_all_options << '\n';
            exit(0);
        }
        po::notify(vm);
    }    
    catch(po::error& e) {
        std::cerr << "ERROR while parsing command line options for " 
                  << filename << ":\n" << e.what() << "\n\n\n"
                  << _M_all_options << '\n';
        exit(1);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR while parsing command line options for " 
                  << filename << ":\n" << e.what() << "\n\n\n"
                  << _M_all_options << '\n';
        exit(1);
    }
}

#endif // __XAVIER_OPTION_PARSE_HPP__
