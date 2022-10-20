#ifndef __XAVIER_STRINGS_HPP__
#define __XAVIER_STRINGS_HPP__

#include <algorithm>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "misc_helper.hpp"

namespace spurt {
    inline std::string number_to_rank(unsigned i) {
        switch (i % 10) {
            case 1:  return std::to_string(i) + "st";
            case 2:  return std::to_string(i) + "nd";
            case 3:  return std::to_string(i) + "rd";
            default: return std::to_string(i) + "th";
        }
    }

    void tokenize(std::vector<std::string>& words,
                  const std::string& phrase,
                  const std::string& delim=" ,.-") {
        char* str=new char[phrase.size()+1];
        std::strcpy(str, phrase.c_str());
        char* pch=strtok (str, delim.c_str());
        while (pch != NULL)
        {
            words.push_back(pch);
            pch = strtok (NULL, delim.c_str());
        }
    }

    std::string& lower_case(std::string& str) {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }

    // Solution by "Eric Malenfant" found on stackoverflow
    unsigned int
    hamming_distance(const std::string& s1, const std::string& s2) {
        if (s1.size()!=s2.size())
            throw std::invalid_argument("hamming distance requires " \
                "string arguments of same length");
        return s1.size() - std::inner_product(s1.begin(), s1.end(),
                                              s2.begin(), 0,
                                              std::plus<>(),
                                              std::equal_to<>());
    }

    // From "Fast, memory efficient Levenshtein algorithm"
    // by Sten Hjelmqvist found referenced on Wikipedia
    // https://en.wikipedia.org/wiki/Levenshtein_distance
    unsigned int
    levenshtein_distance(const std::string& s, const std::string& t) {
        // degenerate cases
        if (s == t) return 0;
        if (s.size() == 0) return t.size();
        if (t.size() == 0) return s.size();

        // create two work vectors of integer distances
        std::vector<int> v0(t.size()+1), v1(t.size()+1);

        // initialize v0 (the previous row of distances)
        // this row is A[0][i]: edit distance for an empty s
        // the distance is just the number of characters to delete from t
        for (unsigned int i=0; i<v0.size(); ++i) v0[i]=i;

        for (unsigned int i=0; i<s.size(); ++i)
        {
            // calculate v1 (current row distances) from the previous row v0

            // first element of v1 is A[i+1][0]
            //   edit distance is delete (i+1) chars from s to match empty t
            v1[0]=i+1;

            // use formula to fill in the rest of the row
            for (unsigned int j=0; j<t.size(); ++j) {
                unsigned int incr=(s[i]==t[j]) ? 0: 1;
                v1[j+1]=minimum<unsigned int>(v1[j]+1, v0[j+1]+1,
                                              v0[j]+incr);
            }

            // copy v1 (current row) to v0 (previous row) for next iteration
            std::swap(v0, v1);
        }
        return v1[t.size()];
    }

    unsigned int
    distance_multistring(std::string arg,
                         const std::vector<std::string>& argv) {
        lower_case(arg);
        std::vector<std::string> lwc(argv.begin(), argv.end());
        std::for_each(lwc.begin(), lwc.end(), [&](std::string& s)
            { lower_case(s);});

        if (argv.size()==1) {
            return arg==lwc[0]? 0: 1;
        }
        else {
            std::string str1(lwc[0]); // words appended
            std::string str2(lwc[0]); // words separated by '_'
            std::string str3(lwc[0]); // words separated by '-'
            std::string str4(lwc[0]); // words separated by ' '
            std::string str5(1, lwc[0][0]); // initials only
            for (size_t i=1; i<lwc.size(); ++i) {
                str1+=    lwc[i];
                str2+='_'+lwc[i];
                str3+='-'+lwc[i];
                str4+=' '+lwc[i];
                str5+=    lwc[i][0];
            }
            std::cout << "arg=" << arg << '\n';
            std::cout << "str1=" << str1 << '\n';
            std::cout << "str2=" << str2 << '\n';
            std::cout << "str3=" << str3 << '\n';
            std::cout << "str4=" << str4 << '\n';
            std::cout << "str5=" << str5 << '\n';
            return minimum(
                levenshtein_distance(arg, str1),
                levenshtein_distance(arg, str2),
                levenshtein_distance(arg, str3),
                levenshtein_distance(arg, str4),
                levenshtein_distance(arg, str5));
        }
    }

    bool match_multistring(std::string arg,
                           const std::vector<std::string>& argv) {
        return distance_multistring(arg, argv)==0;
    }

    bool match_multistring(std::string arg,
                           const std::initializer_list<std::string>&
                           arg_list) {
        std::vector<std::string> argv(arg_list);
        return match_multistring(arg, argv);
    }
} // spurt

#endif
