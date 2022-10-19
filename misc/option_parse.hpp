#ifndef __XAVIER_OPTION_PARSE_HPP__

// STL
#include <exception>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// Boost
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/convert.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

// Add validator-level support for fixed-size arrays to
// boost::program_options
#include <misc/option_parse_array.hpp>
#include <misc/meta_utils.hpp>

namespace boost { namespace program_options {

// Value type to be used when an option requires bounds (lower, upper,
// or both) to be placed on the number of items that a sequence can take.
// This class also covers the case of fixed sized arrays.
template<class T, class charT = char>
class bounded_typed_value : public typed_value<T, charT>
{
public:
    typedef typed_value<T, charT>         base_t;
    typedef bounded_typed_value<T, charT> self_t;

    bounded_typed_value(T* store_to)
        : base_t(store_to), m_min(-1), m_max(-1) {}

    unsigned int min_tokens() const {
        if(m_min < 0) return base_t::min_tokens();
        else return (unsigned int)m_min;
    }

    unsigned int max_tokens() const {
        if(m_max < 0) return base_t::max_tokens();
        else return (unsigned int)m_max;
    }

    self_t* min_tokens(unsigned int min_tokens) {
        if(min_tokens > 1) base_t::multitoken();
        m_min = min_tokens;
        return this;
    }

    self_t* max_tokens(unsigned int max_tokens) {
        if(max_tokens > 1) base_t::multitoken();
        m_max = max_tokens;
        return this;
    }

    self_t* fixed_tokens(unsigned int num_tokens) {
        if(num_tokens > 1) base_t::multitoken();
        m_min = num_tokens;
        m_max = num_tokens;
        return this;
    }

private:
    int m_min;
    int m_max;
};

template<class T, class charT = char>
bounded_typed_value<T, charT>*
bounded_value(T* store_to=0, bool is_required=false)
{
    bounded_typed_value<T, charT>* r=
        new bounded_typed_value<T, charT>(store_to);
    if (is_required) r->required();
    return r;
}


}}

namespace xavier { namespace command_line {

namespace po = boost::program_options;

// Simple class to store general information about a group of options
class option_traits {
public:
    explicit option_traits(bool required=false,
                  bool positional=false,
                  const std::string& group_name="")
        : _default_requirement(required),
          _default_positional(positional),
          _group_name(group_name) {}

    void required(bool required) { _default_requirement = required; }
    void positional(bool positional) { _default_positional = positional; }
    void group_name(const std::string& name) { _group_name = name; }

    bool required() const { return _default_requirement; }
    bool positional() const { return _default_positional; }
    bool has_group() const { return !_group_name.empty(); }
    const std::string& group_name() const { return _group_name; }

protected:
    bool         _default_requirement;
    bool         _default_positional;
    std::string  _group_name;
};

// Helper class to turn a user-defined parsing function
// into a callback function compatible with boost::program_options
// post-parsing notification mechanism
template<typename _Type>
struct custom_parser {
    // for notational convenience
    typedef std::vector<std::string> strvec_t;
    // callback function type expected by boost::program_options
    typedef typename boost::function1<void, strvec_t&> callback_type;
    // user-defined function performing simultaneously validation of a
    // series of strings and value assignment of provided variable
    typedef typename std::function<void (_Type&, const strvec_t&)> parser_type;

    static callback_type callback(_Type& variable,
                                  parser_type parser) {
        return callback_type([&](strvec_t& s) {parser(variable, s);});
    }
};

namespace detail {

// Helper parser that ensures that numbers (and in particular negative
// numbers) won't be parsed as options by removing them from list of
// items to be parsed by other parsers.
std::vector<po::option> skip_numbers(std::vector<std::string>& args)
{
    std::vector<po::option> result;
    int pos = 0;
    while(!args.empty()) {
        const auto& arg = args[0];
        double num;
        if(boost::conversion::try_lexical_convert(arg, num)) {
            result.push_back(po::option());
            po::option& opt = result.back();

            opt.position_key = pos++;
            opt.value.push_back(arg);
            opt.original_tokens.push_back(arg);

            args.erase(args.begin());
        } else {
            break;
        }
    }

    return result;
}

} // detail


/** option_parser: General-purpose parser for command line
    options of arbitrary types and syntax.
    The implementation uses boost::program_options and provides
    a number of useful extensions:
    - simplified interface
    - built-in support for fixed-size array types
    - support for arbitrary types for which the user provides a
      strings to value conversion function.
    - POSIX formatted help message
    - early detection of grammatical errors
 */
class option_parser {
    // notational convenience
    typedef po::options_description          option_group;
    typedef boost::shared_ptr<option_group>  group_ptr;
    typedef std::string                      str_t;

public:
    const str_t default_group_title;

    option_parser(const str_t& program_name,
                  const str_t& synopsis,
                  const str_t& default_section_title="Default Options",
                  size_t line_length=po::options_description::m_default_line_length);

    // Parameters controlling the display of option's arguments
    bool use_default_symbol() const {
        return _use_default_symbol;
    }
    void use_default_symbol(bool use_it=true) {
        _use_default_symbol = use_it;
    }
    bool use_short_symbols() const {
        return _use_short_symbols;
    }
    void use_short_symbols(bool use_it=true) {
        _use_short_symbols = use_it;
        _use_default_symbol = false;
    }
    bool use_brackets() const {
        return _use_brackets;
    }
    void use_brackets(bool use_it=true) {
        _use_brackets = use_it;
    }
    const str_t& default_symbol() const {
        return _default_symbol;
    }
    void default_symbol(const str_t& symbol) {
        _default_symbol = symbol;
    }

    // Add a parameter (default: optional / non-positional)
    template<typename _Type>
    void add_value(const str_t& name, _Type& variable,
                   const str_t& description,
                   const option_traits& traits=option_traits(false),
                   const str_t& symbol="");

    // Add and initialize parameter (default: optional / non-positional)
    template<typename _Type, typename _CompatibleType = _Type>
    void add_value(const str_t& name, _Type& variable,
                   const _CompatibleType& value,
                   const str_t& description,
                   const option_traits& traits=option_traits(false),
                   const str_t& symbol="");

    // Add and initialize boolean flag (ALWAYS: optional / non-positional)
    void add_flag(const str_t& name, bool& variable,
                  const str_t& description,
                  const option_traits& traits);

    // Add fixed-size random access container
    // (e.g., boost::array / std::array)
    // (default: optional / non-positional)
    template<size_t N, typename _Tuple,
             typename _Type = typename _Tuple::value_type>
    void add_tuple(const str_t& name,
                   _Tuple& variable,
                   const str_t& description,
                   const option_traits& traits=option_traits(false),
                   const str_t& symbol="");

    // Add and initialize fixed-size random access container
    // (e.g., boost::array / std::array)
    // By default, the tuple is optional and non-positional
    template<size_t N, typename _Tuple,
             typename _Type = typename _Tuple::value_type,
             typename _CompatibleTuple = _Tuple>
    void add_tuple(const str_t& name,
                   _Tuple& variable,
                   const _CompatibleTuple& value,
                   const str_t& description,
                   const option_traits& traits=option_traits(false),
                   const str_t& symbol="");

    // Add std::vector (default: optional / non-positional)
    template<typename _Type>
    void add_sequence(const str_t& name,
                      std::vector<_Type>& variable,
                      const str_t& description,
                      const option_traits& traits=option_traits(false),
                      const str_t& symbol="");

    // Add std::vector (default: optional / non-positional)
    template<typename _Type, typename _CompatibleType = _Type>
    void add_sequence(const str_t& name,
                      std::vector<_Type>& variable,
                      const std::vector<_CompatibleType>& value,
                      const str_t& description,
                      const option_traits& traits=option_traits(false),
                      const str_t& symbol="");

    // Add std::vector (default: optional / non-positional)
    template<typename _Type>
    void add_bounded_sequence(const str_t& name,
                              std::vector<_Type>& variable,
                              unsigned int min_tok,
                              unsigned int max_tok,
                              const str_t& description,
                              const option_traits& traits=option_traits(false),
                              const str_t& symbol="");

    // Add std::vector (default: optional / non-positional)
    template<typename _Type, typename _CompatibleType = _Type>
    void add_bounded_sequence(const str_t& name,
                              std::vector<_Type>& variable,
                              const std::vector<_CompatibleType>& value,
                              unsigned int min_tok,
                              unsigned int max_tok,
                              const str_t& description,
                              const option_traits& traits=option_traits(false),
                              const str_t& symbol="");

    // Add parameter taking N arguments with custom validation procedure
    // (default: optional / non-positional)
    template<typename _Type, int N>
    void add_custom(const str_t& name,
                    _Type& variable,
                    const str_t& description,
                    typename custom_parser<_Type>::parser_type parser,
                    const option_traits& traits=option_traits(false),
                    const str_t& symbol="");

    // Parse command line options
    void parse(int argc, const char* argv[]);

    // Print selected aspects of the help message
    str_t print_self(bool with_synopsis=true,
                     bool with_usage=true,
                     bool with_options=true,
                     const std::vector<str_t>& with_those_groups =
                     std::vector<str_t>()) const;

private:
    typedef po::positional_options_description positionals_type;

    std::vector<group_ptr>  _groups;
    std::map<str_t, size_t> _group_map;
    positionals_type        _positionals;
    str_t                   _program_name;
    str_t                   _synopsis;
    str_t                   _default_symbol;
    str_t                   _flags_string;
    str_t                   _first_random_length;
    std::vector<str_t>      _options_string;
    std::vector<str_t>      _positionals_string;
    size_t                  _line_length;
    int                     _current_position;
    bool                    _use_default_symbol;
    bool                    _use_brackets;
    bool                    _use_short_symbols;

    // returns usage description with prescribed indentation, whereby
    // each line does not exceed _line_length-1 characters
    str_t usage(size_t initial_indent=0, size_t indent=0) const;

    // returns shortest valid name of an option (with optional dash)
    str_t short_form(const str_t& name, bool with_dash=false);

    // returns symbol representation of option parameter(s)
    str_t untyped_parameter_string(const str_t& s, int N=1) const;
    template<typename _Type>
    str_t typed_parameter_string(const str_t& s, int N=1) const;

    // add boolean switch to usage string
    void add_flag_to_usage_string(const str_t& name);

    // add option description to usage string
    void add_option_to_usage_string(const str_t& name,
                                    const str_t& symbol,
                                    size_t N=1, bool required=false,
                                    bool positional=false);

    // return option group associated with traits
    group_ptr get_group(const option_traits& traits);

    // register positional option
    void register_positional(const str_t& name, int N=1);

    // helper functions to create a po::typed_value object
    template<typename T>
    po::typed_value<T>*
    get_semantic(T& variable, bool required, const str_t& symbol);

    template<typename T>
    po::typed_value<T>*
    get_semantic(T& variable, bool required, const T& value,
                 const str_t& symbol);    //
    //
    // template<>
    // po::typed_value<bool>*
    // get_semantic(bool& variable, bool required, const bool& value,
    //              const str_t& symbol);
    //
    // template<>
    // po::typed_value<char>*
    // get_semantic(char& variable, bool required, const char& value,
    //              const str_t& symbol);
    //
    // template<>
    // po::typed_value<str_t>*
    // get_semantic(str_t& variable, bool required, const str_t& value,
    //              const str_t& symbol);
    //
    // template<>
    // po::typed_value<float>*
    // get_semantic(float& variable, bool required, const float& value,
    //              const str_t& symbol);
    //
    // template<>
    // po::typed_value<double>*
    // get_semantic(double& variable, bool required, const double& value,
    //              const str_t& symbol);

    template<typename T>
    po::typed_value<T>*
    get_semantic(T& variable, bool required, const T& value,
                 const str_t& value_as_string,
                 const str_t& symbol);

    // reflow a paragraph and apply prescribed indentation
    str_t reflow(const str_t& s, size_t init_indent=0,
                       size_t indent=0) const;
};

template<typename T, typename Enable=void>
struct string_converter {
    static std::string convert(const T& value) {

        // std::cout << "called vanilla string_converter for " << value << '\n';

        return boost::lexical_cast<std::string>(value);
    }
};

template<typename T>
struct string_converter<T, typename std::enable_if<std::is_floating_point<T>::value>::type> {
    static std::string convert(const T& value) {

        // std::cout << "string converter called for floating point value: "
        //     << value << std::endl;

        std::ostringstream os;
        os << /* std::fixed << */ std::setprecision(2) << value;
        return os.str();
    }
};


} // command_line
} // xavier

std::ostream& operator<<(std::ostream& oss,
                         const xavier::command_line::option_parser& parser) {
    oss << parser.print_self();
    return oss;
}

inline xavier::command_line::option_parser::
option_parser(const str_t& program_name,
              const str_t& synopsis,
              const str_t& default_section_title,
              size_t line_length)
    : default_group_title(default_section_title),
      _groups(), _positionals(), _synopsis(synopsis),
      _default_symbol("arg"),
      _line_length(line_length),
      _current_position(0),
      _use_default_symbol(false), _use_brackets(true),
      _use_short_symbols(true) {

    boost::filesystem::path p(program_name);
    _program_name = p.filename().string();

    // create a group for general options

    option_group* general =
        new option_group(default_group_title);
    _groups.push_back(boost::shared_ptr<option_group>(general));
    _groups.back()->add_options()("help,h", "Print this message");
    _group_map[default_group_title] = 0;
    _flags_string += "-h";
}

inline std::string xavier::command_line::option_parser::
print_self(bool with_synopsis, bool with_usage, bool with_options,
           const std::vector<str_t>& with_those_groups) const {
    std::ostringstream oss;
    if (with_synopsis)
        oss << reflow(_program_name + ": " + _synopsis) << "\n\n";
    if (with_usage)
        oss << usage() << "\n\n";
    if (with_options) {
        option_group selected_groups(_line_length);
        if (with_those_groups.empty()) {
            for (auto it = _groups.begin(); it!=_groups.end() ;
                 ++it) {
                selected_groups.add(**it);
            }
        }
        else {
            for (auto it1 = with_those_groups.begin() ;
                 it1!=with_those_groups.end() ; ++it1) {
                auto it2 = _group_map.find(*it1);
                if (it2 != _group_map.end())
                    selected_groups.add(*_groups[it2->second]);
            }
        }
        oss << selected_groups;
    }
    return oss.str();
}

template<typename T>
inline boost::program_options::typed_value<T>*
xavier::command_line::option_parser::
get_semantic(T& variable, bool required, const str_t& symbol) {
    po::typed_value<T>* v =
        po::value<T>(&variable)->value_name(symbol.c_str());
    if (required) v->required();
    return v;
}

template<typename T>
inline boost::program_options::typed_value<T>*
xavier::command_line::option_parser::
get_semantic(T& variable, bool required, const T& value,
             const str_t& symbol) {
    po::typed_value<T>* v = get_semantic(variable, required, symbol);
    return v->default_value(value);
}

namespace xavier { namespace command_line {
template<>
inline boost::program_options::typed_value<bool>*
option_parser::
get_semantic(bool& variable, bool required, const bool& value,
             const str_t& symbol) {
    str_t value_as_string = (value ? "true" : "false");
    return get_semantic<bool>(variable, required, value, value_as_string,
                              symbol);
}

template<>
inline boost::program_options::typed_value<char>*
option_parser::
get_semantic(char& variable, bool required, const char& value,
             const str_t& symbol) {
    str_t value_as_string = "\'" + str_t(1, value) + "\'";
    return get_semantic<char>(variable, required, value, value_as_string,
                              symbol);
}

template<>
inline boost::program_options::typed_value<
    xavier::command_line::option_parser::str_t>*
option_parser::
get_semantic(str_t& variable, bool required, const str_t& value,
             const str_t& symbol) {
    str_t value_as_string = "\"" + value + "\"";
    return get_semantic<str_t>(variable, required, value, value_as_string,
                               symbol);
}

template<>
inline boost::program_options::typed_value<float>*
option_parser::
get_semantic(float& variable, bool required, const float& value,
             const str_t& symbol) {
    std::ostringstream os;
    os << std::setprecision(2);
    if (fabs(value)<0.01) os << std::scientific << value;
    else os << std::fixed << value;
    return get_semantic<float>(variable, required, value, os.str(),
                               symbol);
}

template<>
inline boost::program_options::typed_value<double>*
option_parser::
get_semantic(double& variable, bool required, const double& value,
             const str_t& symbol) {
    std::ostringstream os;
    os << std::setprecision(2);
    if (fabs(value)<0.01) os << std::scientific << value;
    else os << std::fixed << value;
    return get_semantic<double>(variable, required, value, os.str(),
                                symbol);
}
} // command_line
} // xavier

template<typename T>
inline boost::program_options::typed_value<T>*
xavier::command_line::option_parser::
get_semantic(T& variable, bool required, const T& value,
             const str_t& value_as_string,
             const str_t& symbol) {
    po::typed_value<T>* v = get_semantic(variable, required, symbol);
    return v->default_value(value, value_as_string.c_str());
}

inline void xavier::command_line::option_parser::
register_positional(const str_t& name, int N) {
    if (_current_position == -1) {
        std::ostringstream os;
        os << "Invalid sequence of arguments: "
           << "Attempting to add <" << name << "> as positional argument after "
           << "arbitrary length "
           << "parameter <" << _first_random_length << ">.";
        throw std::runtime_error(reflow(os.str()));
    }
    if (N == -1) {
        _first_random_length = name;
        _current_position = -1;
    }
    else _current_position += N;
    _positionals.add(name.c_str(), _current_position);
}

inline std::string xavier::command_line::option_parser::
reflow(const str_t& str, size_t init_indent, size_t indent) const {
    typedef boost::char_separator<char>    separator_t;
    typedef boost::tokenizer<separator_t>  tokenizer_t;

    separator_t sep(" ");
    tokenizer_t tokens(str, sep);

    str_t indent_1, indent_2;
    indent_1.assign(init_indent, ' ');
    indent_2.assign(init_indent + indent, ' ');

    std::ostringstream oss;
    oss << indent_1;
    size_t line_length = init_indent;
    for (tokenizer_t::iterator it=tokens.begin() ; it!=tokens.end(); ++it) {
        if (it == tokens.begin()) {
            oss << *it;
            line_length += it->size();
        }
        else if (it->size() + line_length + 1 < _line_length) {
            oss << " " << *it;
            line_length += it->size() + 1;
        }
        else {
            oss << '\n' << indent_2 << *it;
            line_length = it->size() + indent_2.size();
        }
    }

    return oss.str();
}

// Display grammar for option <name> with value symbol <sym>
// where <short> =
//    * "-<c>" if <c> is the one-character short name of the option
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
usage(size_t initial_indent, size_t indent) const {

    str_t indent_1, indent_2;
    indent_1.assign(initial_indent, ' ');
    indent_2.assign(initial_indent + indent, ' ');

    std::ostringstream oss;
    oss << indent_1 << _program_name << " [" << _flags_string << "]";
    size_t length = initial_indent + _program_name.size() +
                    _flags_string.size() + 1;
    for (size_t i=0 ; i<_options_string.size() ; ++i) {
        const str_t& item = _options_string[i];
        if (item.size() + length + 1 < _line_length) {
            oss << " " << item;
            length += item.size() + 1;
        }
        else {
            oss << '\n' << indent_2 << item;
            length = item.size() + indent + initial_indent;
        }
    }
    for (size_t i=0 ; i<_positionals_string.size() ; ++i) {
        const str_t& item = _positionals_string[i];
        if (item.size() + length + 1 < _line_length) {
            oss << " " << item;
            length += item.size() + 1;
        }
        else {
            oss << '\n' << indent_2 << item;
            length = item.size() + initial_indent + indent;
        }
    }
    return oss.str();
}

// <short>
inline std::string xavier::command_line::option_parser::
short_form(const str_t& name, bool with_dash) {
    size_t n = name.find(',');
    str_t s;
    if (n != str_t::npos) {
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
untyped_parameter_string(const str_t& s, int N) const {
    if (N < 0) {
        return s + "...";
    }
    else if (N == 1) return s;
    else if (N*s.size() + N - 1 <= 10) { // somewhat random
        str_t r = s;
        for (size_t i=1 ; i<N; ++i) {
            r += " " + s;
        }
        return r;
    }
    else {
        return boost::lexical_cast<str_t>(N) + " " + s + "'s";
    }
}

template<typename _Type>
inline std::string xavier::command_line::option_parser::
typed_parameter_string(const str_t& s, int N) const {
    if (!s.empty()) return s;
    else if (_use_default_symbol)
        return untyped_parameter_string(_default_symbol, N);
    else {
        return untyped_parameter_string(type_as_string<_Type>
                                            (_use_brackets,
                                             _use_short_symbols),
                                        N);
    }
}

inline void xavier::command_line::option_parser::
add_flag_to_usage_string(const str_t& name) {
    str_t s = short_form(name, false);
    if (s.size() == 1) _flags_string += s;
    else _options_string.push_back("[" + short_form(name, true) + "]");
}

inline void xavier::command_line::option_parser::
add_option_to_usage_string(const str_t& name,
                           const str_t& val_str,
                           size_t N, bool required, bool positional) {
    str_t opt_str = short_form(name, true) + " " + val_str;
    str_t form_str;
    if (positional) form_str = val_str + "|" + opt_str;
    else form_str = opt_str;
    str_t usage_str;
    if (required) usage_str = form_str;
    else usage_str = "[" + form_str + "]";
    if (positional) _positionals_string.push_back(usage_str);
    else _options_string.push_back(usage_str + " ");
}

inline xavier::command_line::option_parser::group_ptr
xavier::command_line::option_parser::
get_group(const option_traits& traits) {
    const str_t& name = traits.group_name();
    if (!name.empty()) {
        std::map<str_t, size_t>::const_iterator it;
        it = _group_map.find(name);
        if (it != _group_map.end()) {
            // fetch existing group
            return _groups[it->second];
        }
        else {
            // create new group
            option_group* new_group =
                new option_group(name);
            _groups.push_back(group_ptr(new_group));
            _group_map[name] = _groups.size()-1;
            return _groups.back();
        }
    }
    else { // default group
        return _groups[0];
    }
}

template<typename _Type>
inline void xavier::command_line::option_parser::
add_value(const str_t& name, _Type& variable,
          const str_t& description,
          const option_traits& traits,
          const str_t& symbol) {

    const str_t s = typed_parameter_string<_Type>(symbol);
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
add_value(const str_t& name, _Type& variable,
          const _CompatibleType& value,
          const str_t& description,
          const option_traits& traits,
          const str_t& symbol) {

    _Type init_value = static_cast<_Type>(value);
    const str_t s = typed_parameter_string<_Type>(symbol);
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
add_flag(const str_t& name, bool& variable,
         const str_t& description,
         const option_traits& traits) {

    variable = false;
    get_group(traits)->add_options()(
        name.c_str(), po::bool_switch(&variable), description.c_str()
    );
    add_flag_to_usage_string(name);
}

template<size_t N, typename _Tuple, typename _Type, typename _CompatibleTuple>
inline void xavier::command_line::option_parser::
add_tuple(const str_t& name, _Tuple& variable,
          const _CompatibleTuple& value,
          const str_t& description,
          const option_traits& traits,
          const str_t& symbol) {

    _Tuple initial_value;
    // std::ostringstream os1, os2;
    // os1 << "(";
    // os2 << "(";
    for (size_t i=0 ; i<N ; ++i) {
        // os1 << value[i] << (i<N-1 ? "," : ")");
        initial_value[i] = static_cast<_Type>(value[i]);
        // os2 << initial_value[i] << (i<N-1 ? "," : ")");
    }
    //
    // std::cout << "input value=" << os1.str()
    //     << ", initial_value=" << os2.str() << '\n';

    // std::ostringstream os;
    // os << initial_value;
    // str_t value_as_string = os.str();
    str_t value_as_string = "(";
    for (size_t i=0 ; i<N ; ++i) {
        value_as_string += string_converter<_Type>::convert(value[i]);
        // std::cout << "added " << value[i] << " to " << value_as_string << '\n';
        if (i < N-1) value_as_string += ",";
        else value_as_string += ")";
    }

    const str_t s = typed_parameter_string<_Type>(symbol, N);

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
add_tuple(const str_t& name, _Tuple& variable,
          const str_t& description,
          const option_traits& traits,
          const str_t& symbol) {

    // std::cout << "entering add_tuple #1\n";

    const str_t s = typed_parameter_string<_Type>(symbol, N);
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
add_sequence(const str_t& name, std::vector<_Type>& variable,
             const str_t& description,
             const option_traits& traits,
             const str_t& symbol) {

    const str_t s = typed_parameter_string<_Type>(symbol, -1);
    get_group(traits)->add_options()(
        name.c_str(),
        get_semantic(variable, traits.required(), s)->multitoken(),
        description.c_str()
    );
    add_option_to_usage_string(name, s, -1, traits.required(),
                               traits.positional());
    if (traits.positional()) register_positional(name, -1);
    else {
        if (_current_position >= 0) _first_random_length = name;
        _current_position = -1; // no positional possible after a sequence
    }
}

template<typename _Type, typename _CompatibleType>
inline void xavier::command_line::option_parser::
add_sequence(const str_t& name,
             std::vector<_Type>& variable,
             const std::vector<_CompatibleType>& value,
             const str_t& description,
             const option_traits& traits,
             const str_t& symbol) {
    std::vector<_Type> init_value(value.size());
    for (size_t i=0; i<value.size(); ++i) {
        init_value[i] = static_cast<_Type>(value[i]);
    }

    str_t value_as_string = "(";
    for (size_t i=0 ; i<init_value.size(); ++i) {
        value_as_string += string_converter<_Type>::convert(init_value[i]);
        // std::cout << "added " << value[i] << " to " << value_as_string << '\n';
        if (i < init_value.size()-1) value_as_string += ",";
        else value_as_string += ")";
    }

    const str_t s = typed_parameter_string<_Type>(symbol, -1);
    get_group(traits)->add_options()(
        name.c_str(),
        get_semantic(variable, traits.required(), init_value,
                     value_as_string, s)->multitoken(),
        description.c_str()
    );
    add_option_to_usage_string(name, s, -1, traits.required(),
                               traits.positional());
    if (traits.positional()) register_positional(name, -1);
    else {
        if (_current_position >= 0) _first_random_length = name;
        _current_position = -1; // no positional possible after a sequence
    }
}

template<typename _Type>
inline void xavier::command_line::option_parser::
add_bounded_sequence(const str_t& name,
                     std::vector<_Type>& variable,
                     unsigned int min_tok,
                     unsigned int max_tok,
                     const str_t& description,
                     const option_traits& traits,
                     const str_t& symbol) {

    const str_t s = typed_parameter_string<_Type>(symbol, -1);
    get_group(traits)->add_options()(
        name.c_str(),
        po::bounded_value< std::vector<_Type> >(&variable, traits.required())
            ->min_tokens(min_tok)
            ->max_tokens(max_tok)
            ->value_name(s.c_str()),
        description.c_str()
    );
    add_option_to_usage_string(name, s, -1, traits.required(),
                               traits.positional());
    if (traits.positional()) register_positional(name, -1);
    else {
        if (_current_position >= 0) _first_random_length = name;
        _current_position = -1; // no positional possible after a sequence
    }
}

template<typename _Type, typename _CompatibleType>
inline void xavier::command_line::option_parser::
add_bounded_sequence(const str_t& name, std::vector<_Type>& variable,
                     const std::vector<_CompatibleType>& value,
                     unsigned int min_tok,
                     unsigned int max_tok,
                     const str_t& description,
                     const option_traits& traits,
                     const str_t& symbol) {
    std::vector<_Type> init_value(value.size());
    for (size_t i=0; i<value.size(); ++i) {
        init_value[i] = static_cast<_Type>(value[i]);
    }

    str_t value_as_string = "(";
    for (size_t i=0 ; i<init_value.size(); ++i) {
        value_as_string += string_converter<_Type>::convert(init_value[i]);
        // std::cout << "added " << value[i] << " to " << value_as_string << '\n';
        if (i < init_value.size()-1) value_as_string += ",";
        else value_as_string += ")";
    }

    const str_t s = typed_parameter_string<_Type>(symbol, -1);
    get_group(traits)->add_options()(
        name.c_str(),
        po::bounded_value< std::vector<_Type> >(&variable, traits.required())
            ->min_tokens(min_tok)
            ->max_tokens(max_tok)
            ->value_name(s.c_str())
            ->default_value(init_value, value_as_string.c_str()),
        description.c_str()
    );
    add_option_to_usage_string(name, s, -1, traits.required(),
                               traits.positional());
    if (traits.positional()) register_positional(name, -1);
    else {
        if (_current_position >= 0) _first_random_length = name;
        _current_position = -1; // no positional possible after a sequence
    }
}

template<typename _Type, int N>
inline void xavier::command_line::option_parser::
add_custom(const str_t& name, _Type& variable,
           const str_t& description,
           typename custom_parser<_Type>::parser_type parser,
           const option_traits& traits,
           const str_t& symbol) {

    typedef custom_parser<_Type> parser_type;

    const std::string s = typed_parameter_string<_Type>(symbol);
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
option_parser::parse(int argc, const char* argv[]) {
    // consolidate all option groups into a single overall group
    option_group all_options(_line_length);
    for (auto it = _groups.begin(); it!=_groups.end() ;
         ++it) {
        all_options.add(**it);
    }

    po::variables_map vm;
    try {
        po::store(
            po::command_line_parser(argc, argv)
            .extra_style_parser(&detail::skip_numbers)
            .options(all_options)
            .positional(_positionals).run(),
            vm
        );

        if (vm.count("help")) {
            std::cout << print_self() << '\n';
            exit(1);
        }
        po::notify(vm);
    }
    catch(po::error& e) {
        std::cerr << "ERROR while parsing command line options for "
                  << _program_name << ":\n" << e.what() << "\n\n\n"
                  << print_self() << '\n';
        exit(0);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR while parsing command line options for "
                  << _program_name << ":\n" << e.what() << "\n\n\n"
                  << print_self() << '\n';
        exit(0);
    }
}

#endif // __XAVIER_OPTION_PARSE_HPP__
