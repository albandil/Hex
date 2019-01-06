//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_DB_UTIL_H
#define HEX_DB_UTIL_H

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "hex-arrays.h"
#include "hex-vec3d.h"

#include "db.h"

/**
 * @brief Energy units change.
 * 
 * Returns factor that can be used to transform from the unit system A
 * to the unit system B.
 */
double change_units (eUnit A, eUnit B);

/**
 * @brief Lengths units change.
 * 
 * Returns factor that can be used to transform from the unit system A
 * to the unit system B.
 */
double change_units (lUnit A, lUnit B);

/**
 * @brief Angular units change.
 * 
 * Returns factor that can be used to transform from the unit system A
 * to the unit system B.
 */
double change_units (aUnit A, aUnit B);

/**
 * @brief Energy unit name.
 * 
 * Return energy unit name as string.
 */
std::string unit_name (eUnit u);

/**
 * @brief Length unit name.
 * 
 * Return length unit name as string.
 */
std::string unit_name (lUnit u);

/**
 * @brief Length unit name.
 * 
 * Return length unit name as string.
 */
std::string unit_name (aUnit u);

/**
 * Write out std::pair.
 */
inline std::ostream & operator << (std::ostream & os, std::pair<geom::vec3d,geom::vec3d> const & p)
{
    os << p.first << " " << p.second;
    return os;
}

/**
 * Read in std::pair.
 */
inline std::istream & operator >> (std::istream & is, std::pair<geom::vec3d,geom::vec3d> & p)
{
    is >> p.first;
    is >> p.second;
    return is;
}

/**
 * Read in std::vector&lt;int&gt;.
 * Allowed are comma-separated values and integer ranges.
 */
inline std::istream & operator >> (std::istream & is, iArray & p)
{
    std::string token;
    p.resize(0);
    while (std::getline(is, token, ','))
    {
        std::size_t idx, idx1, idx2;
        int value = std::stoi(token, &idx);

        // convert token to integer
        if (idx == token.size())
        {
            p.push_back(value);
        }

        // convert token to range
        else
        {
            std::size_t pos = token.find('-', 1);
            if (pos == std::string::npos)
                HexException("Failed to parse \"%s\"", token.c_str());

            std::string str_a = token.substr(0, pos);
            std::string str_b = token.substr(pos + 1);
            int a = std::stoi(str_a, &idx1);
            int b = std::stoi(str_b, &idx2);
            if (idx1 != str_a.size() or idx2 != str_b.size())
                HexException("Failed to parse \"%s\"", token.c_str());

            iArray values = linspace(a, b, b - a + 1);
            p.append(values.begin(), values.end());
        }
    }
    return is;
}

/**
 * Read data from standard input.
 */
template<typename T> std::vector<T> readStandardInput ()
{
    std::vector<T> data;

    T x;
    while (not std::cin.eof())
    {
        std::cin >> std::ws;
        std::cin >> x;
        std::cin >> std::ws;
        data.push_back(x);
    }

    return data;
}

/**
 * @brief Convert dictionary entry to a numeric type.
 * 
 * Being given a dictionary (= string-string map) and a keyword,
 * the function finds a correct entry and returns its value converted
 * to the template datatype.
 * @param dict Dictionary to search in.
 * @param keyword Entry to look for.
 * @param name Identification of the calling authority for use in error
 *             message if the entry is not find.
 */
template <typename T> T Conv
(
    std::map<std::string,std::string> const & dict,
    std::string const & keyword,
    std::string const & name
)
{
    // check existence of the keyword
    std::map<std::string,std::string>::const_iterator it = dict.find(keyword);
    if (it == dict.end())
        throw exception ("ERROR: \"%s\" requires specifying the parameter \"--%s\"!\n", name.c_str(), keyword.c_str());

    // convert to int
    T x;
    std::istringstream ss(it->second);
    ss >> x;
    return x;
}

#endif
