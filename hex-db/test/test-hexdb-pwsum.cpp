//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2026, Jakub Benda, Charles University in Prague                    //
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

// Test of the summation of the partial ionization cross sections over the total
// angular momentum, which 'ics' performs when asked for a negative partial wave.
//
// The partial cross sections decay geometrically, so what lies past the last computed
// one is added as the remainder of a geometric progression rather than dropped. The
// test builds a family of radial amplitudes whose partial cross sections are exactly
// geometric,
//
//     Xi_L(u) = r^(L/2) G(u)   =>   sigma_L = sigma_0 r^L ,
//
// for which the extrapolation is not merely a good estimate but exact, whatever
// partial wave it is applied at:
//
//     sum(0..n) sigma_L + sigma_n q/(1-q)
//         = sigma_0 (1 - r^(n+1))/(1 - r) + sigma_0 r^(n+1)/(1 - r)
//         = sigma_0 / (1 - r) .
//
// So the summed cross section has to come out as sigma_0/(1-r) to the precision the
// table is printed with, and the checks below say so. sigma_0 itself is taken from
// the same code, by asking for the partial wave zero; what is under test here is the
// summation, and the absolute normalization of a single partial wave is pinned by
// test-hexdb-ionization instead.
//
// The second family is there for the other half of the guard. Ionization transfers at
// least the ionization energy, so the adiabatic cutoff bounds the ratio by
// exp(-Eion/Ei); a sequence decaying more slowly than that is not the geometric tail
// of a converging partial wave series but the noise floor of the amplitudes, where the
// factor q/(1-q) would multiply round-off by an arbitrarily large number. The
// extrapolation has to refuse it and say so.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "hex-arrays.h"
#include "hex-chebyshev.h"

#include "../src/db.h"
#include "../src/quantities.h"

// --------------------------------------------------------------------------------- //

const int testNi = 1;
const int testLi = 0;
const int testMi = 0;

// A family of partial cross sections, all stored at one energy so that a single
// database can hold several of them.
//
//   - 'convergent' decays fast enough for the adiabatic bound exp(-Eion/Ei) to accept
//     it and keeps the ratio up to the last partial wave stored;
//   - 'too slow' decays more slowly than the bound allows. That is not the tail of a
//     converging series but the noise floor of the amplitudes, and extrapolating it
//     would multiply round-off by q/(1-q); it has to be refused outright;
//   - 'noise floor' is the realistic mixture: it decays like the first up to Lgeom and
//     then stops decaying. The extrapolation has to break off at Lgeom and discard the
//     partial waves beyond it rather than add them in.
struct Family
{
    double Ei;              // impact energy [Ry]
    double ratio;           // ratio of two successive partial cross sections
    int Lgeom;              // highest partial wave that still follows that ratio
    int Lmax;               // highest partial wave stored
    double bound;           // exp(-Eion/Ei), quoted in the output
    char const * name;
};

const Family families [] =
{
    { 5.0, 0.50, 6,  6, 0.8187, "convergent"  },
    { 8.0, 0.95, 6,  6, 0.8825, "too slow"    },
    { 6.0, 0.50, 6, 10, 0.8465, "noise floor" }
};

const double tolerance = 1e-5;

// --------------------------------------------------------------------------------- //

// Radial amplitude of the partial wave zero, as a function of the normalized momentum
// u = k1/kmax. Its shape does not matter, only that it is not degenerate.
Complex G (double u)
{
    return Complex(1.0, 0.4) * (1.0 + 0.5 * u * u);
}

std::string str (double x)
{
    std::ostringstream oss;
    oss.precision(17);
    oss << x;
    return oss.str();
}

// Run a scattering quantity and return its output table as text.
std::string query (std::string const & quantity, std::map<std::string,std::string> const & params)
{
    std::ostringstream out;
    std::istringstream in ("");

    std::streambuf * cout_buf = std::cout.rdbuf(out.rdbuf());
    std::streambuf * cin_buf  = std::cin.rdbuf(in.rdbuf());

    try
    {
        hex_run({ quantity }, params);
    }
    catch (std::exception const & e)
    {
        std::cout.rdbuf(cout_buf);
        std::cin.rdbuf(cin_buf);
        std::cerr << "Query of \"" << quantity << "\" failed: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::cout.rdbuf(cout_buf);
    std::cin.rdbuf(cin_buf);

    return out.str();
}

// The fields of the single data row of a table.
std::vector<std::string> row (std::string const & table)
{
    std::vector<std::vector<std::string>> data;

    std::istringstream iss (table);
    std::string line;

    while (std::getline(iss, line))
    {
        std::istringstream fields (line);
        std::vector<std::string> fs;
        std::string field;

        while (fields >> field)
            fs.push_back(field);

        if (not fs.empty() and fs.front().front() != '#')
            data.push_back(fs);
    }

    if (data.size() != 1)
    {
        std::cerr << "Expected one data row, got:\n" << table << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return data.front();
}

// Ask 'ics' for one partial wave, or for their sum when 'ell' is negative. Returns
// the cross section; the estimate of the remainder goes to 'tail'.
double ics (double Ei, int S, int ell, double * tail = nullptr)
{
    std::map<std::string,std::string> p =
    {
        { "ni", std::to_string(testNi) }, { "li", std::to_string(testLi) },
        { "mi", std::to_string(testMi) }, { "nf", "0" }, { "lf", "0" },
        { "mf", "0" }, { "S", std::to_string(S) },
        { "ell", std::to_string(ell) }, { "Ei", str(Ei) }
    };

    std::vector<std::string> fields = row(query("ics", p));

    if (tail != nullptr)
        *tail = std::stod(fields.at(2));

    return std::stod(fields.at(1));
}

// --------------------------------------------------------------------------------- //

int failures = 0;

void check (std::string const & what, double found, double expected)
{
    double error = std::abs(found - expected) / std::max(std::abs(expected), 1e-300);
    bool ok = (error < tolerance);

    std::cout << (ok ? "  ok   " : "  FAIL ") << what
              << ": got " << found << ", expected " << expected
              << " (relative error " << error << ")" << std::endl;

    if (not ok)
        failures++;
}

void require (std::string const & what, bool condition)
{
    std::cout << (condition ? "  ok   " : "  FAIL ") << what << std::endl;

    if (not condition)
        failures++;
}

// --------------------------------------------------------------------------------- //

int main (void)
{
    Eunits = eUnit_Ry;
    Lunits = lUnit_au;
    Aunits = aUnit_rad;

    hex_initialize("hex-pwsum.db");
    hex_new();

    // Write the two families the way hex-ecs writes its expansions. One angular
    // channel per total angular momentum is enough; the summation does not look
    // inside a partial wave.
    std::ofstream sql ("batch-pwsum.sql");
    sql << "BEGIN TRANSACTION;\n";

    for (Family const & f : families)
    {
        double kmax = std::sqrt(f.Ei - 1.0/(testNi*testNi));

        for (int L = 0; L <= f.Lmax; L++)
        {
            // sigma_L is quadratic in the amplitude, hence the square root; past
            // Lgeom the scale is held fixed, which is a flat sequence of ratio one
            double scale = std::pow(f.ratio, 0.5 * std::min(L, f.Lgeom));

            Chebyshev<double,Complex> cb;
            cb.generate([&](double k1) { return scale * G(k1 / kmax); }, 64, 0., kmax);

            // one angular channel per total angular momentum is enough; the
            // summation does not look inside a partial wave
            for (int S = 0; S <= 1; S++)
            {
                sql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                    << testNi << "," << testLi << "," << testMi << ","
                    << L << "," << S << "," << f.Ei << ","
                    << L << ",0," << cb.coeffs().toBlob().c_str() << ");\n";
            }
        }
    }

    sql << "COMMIT;\n";
    sql.close();

    hex_import("batch-pwsum.sql");
    hex_update();

    // ------------------------------------------------------------------------- //

    for (Family const & f : families)
    {
        std::cout << f.name << ": Ei = " << f.Ei << " Ry, ratio " << f.ratio
                  << ", partial waves 0 to " << f.Lmax
                  << " (geometric to " << f.Lgeom << ")"
                  << ", adiabatic bound " << f.bound << std::endl;

        bool refused = (f.ratio > f.bound);

        for (int S = 0; S <= 1; S++)
        {
            std::string tag = ", S = " + std::to_string(S);

            double sigma0 = ics(f.Ei, S, 0);
            double tail = 0;
            double total = ics(f.Ei, S, -1, &tail);

            // partial sum of the progression through the partial wave n
            auto partial = [&](int n) -> double
            {
                return sigma0 * (1.0 - std::pow(f.ratio, n + 1)) / (1.0 - f.ratio);
            };

            if (refused)
            {
                // nothing may be extrapolated: the bare sum of everything stored
                check("total" + tag, total, partial(f.Lmax));
                require("no remainder reported" + tag, std::isnan(tail));
                continue;
            }

            // the extrapolation is exact for a geometric progression, wherever it is
            // applied, so the whole series must come out however many terms are stored
            check("total" + tag, total, sigma0 / (1.0 - f.ratio));

            require("a remainder was added" + tag, std::isfinite(tail) and tail > 0.);

            // and it must have been applied at the last partial wave that still
            // follows the ratio -- not earlier, which would waste the data that is
            // there, and not later, which would drag the noise floor into the sum
            check("truncation point" + tag, total - tail, partial(f.Lgeom));
        }
    }

    // ------------------------------------------------------------------------- //

    std::cout << (failures == 0 ? "All tests passed." : "Some tests failed.") << std::endl;

    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
