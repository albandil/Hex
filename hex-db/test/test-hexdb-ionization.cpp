//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

// Regression test of the ionization quantities derived from the 'ionf' table.
//
// The test fills a database with synthetic radial ionization amplitudes Ξ, whose
// Chebyshev expansions are written exactly the way hex-ecs writes them, and then
// checks that every quantity computed from them stays consistent with
//
//                                       π/4
//                    LS      2S + 1      ⌠         L         2
//                   σ     = --------  ⎮   ∑  |Ξ      (kₘₐₓ sin β)|  dβ  ,
//                    ion     4 kᵢ       ⌡  ℓ₁ℓ₂   ℓ₁ℓ₂
//                                       0
//
// which is what the 'ics' quantity evaluates and what an independent derivation of
// the hyperspherical surface-integral extraction in hex-ecs gives. The checks are
// written so that each one fails on its own:
//
//   1. 'ics' against that formula, for both spins and all three total angular momenta;
//   2. 'sdcsE' against dσ/dx = σ' · E_tot/(2 k₁k₂), and 'ionf' against Ξ/(k₁k₂),
//      both pointwise;
//   3. ∫ 'sdcsE' dx = ∑_L 'ics', an invariant between two quantities;
//   4. ∫ 'sdcs12' d(cos θ₁₂) = ∑_L 'ics', likewise;
//   5. 'sdcs12' against pinned values, which the other checks cannot see because
//      they only probe λ = 0, where the relative Coulomb phases cancel;
//   6. 'tdcs' = (2S+1)/4 |'ionamp'|² k₁k₂/√Eᵢ, tying the two amplitude-valued
//      quantities to the cross section built from the same expansions.
//
// The synthetic Ξ respect the exchange symmetry Ξ_{ℓ₁ℓ₂}(k₁) = ±Ξ_{ℓ₂ℓ₁}(k₂) of a
// real solution, which checks 4 and 5 rely on: 'sdcs12' integrates over the k₁ > k₂
// half of the hyperangle and 'ics' over the other one.
//
// Not covered: the hex-ecs side of the same convention, that is, the spin weight in
// Amplitudes::computeSigmaIon_ and the expansions that Amplitudes::writeSQL_files
// emits. Reaching either needs a solved wave function, so check 1 pins the shared
// convention here and hex-ecs has to be held to it by its own tests.

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

// --------------------------------------------------------------------------------- //

// the scattering event the whole test is built around
const int    testNi = 1;
const int    testLi = 0;
const int    testMi = 0;
const double testEi = 5.0;                     // impact energy [Ry]
const double testEtot = testEi - 1.0/(testNi*testNi);   // energy left to the two electrons
const double testKmax = std::sqrt(testEtot);   // 2, deliberately not 1, so that an
                                               // expansion misread as living on [0,1]
                                               // instead of [0,kₘₐₓ] shows up
const double testKi = std::sqrt(testEi);

// Relative tolerance of the comparisons. The quantities print through OutputTable,
// which uses the default six significant digits, so a parsed value cannot be trusted
// much beyond 1e-6; the quadratures below add a little on top of that. This is still
// three orders of magnitude tighter than the smallest normalization error the test is
// meant to catch, a factor of 4/3.
const double tolerance = 1e-4;

// --------------------------------------------------------------------------------- //

// One angular channel of the synthetic solution. The radial amplitude is given as a
// function of the normalized momentum u = k₁/kₘₐₓ ∈ [0,1], which makes the exchange
// symmetry easy to state: the partner channel (ℓ₂,ℓ₁) carries G(√(1-u²)).
struct Channel
{
    int L, l1, l2;
    Complex (*G) (double u);
};

// diagonal channel, symmetric under u ↔ √(1-u²) because it depends on u only
// through (u² - ½), which changes sign under the exchange
Complex G00 (double u)
{
    double t = u*u - 0.5;
    return Complex(1.0, 0.4) * (1.0 + 2.0*t*t);
}

// off-diagonal channel, arbitrary
Complex G01 (double u)
{
    return Complex(0.7, -0.3) * (1.0 + u);
}

// its exchange partner
Complex G10 (double u)
{
    return G01(std::sqrt(std::max(0.0, 1.0 - u*u)));
}

// a second diagonal channel, symmetric for the same reason as G00
Complex G11 (double u)
{
    double t = u*u - 0.5;
    return Complex(0.6, -0.5) * (1.0 + 0.5*t*t);
}

// The (1,1) channel is what makes the bipolar harmonic non-trivial: with ℓ₁ = 0 or
// ℓ₂ = 0 the Clebsch-Gordan coefficient leaves only m = 0 in the sum, and a wrong
// magnetic quantum number in the second spherical harmonic cannot be seen.
const std::vector<Channel> channels =
{
    { 0, 0, 0, &G00 },
    { 1, 0, 1, &G01 },
    { 1, 1, 0, &G10 },
    { 2, 1, 1, &G11 }
};

// --------------------------------------------------------------------------------- //

// Format a parameter for the quantity. Not std::to_string, which writes six decimal
// places and would both coarsen the quadrature nodes below and round the smallest
// energy share down to zero.
std::string str (double x)
{
    std::ostringstream oss;
    oss.precision(17);
    oss << x;
    return oss.str();
}

std::string str (int i)
{
    return std::to_string(i);
}

// Run a scattering quantity and return its output table as text. The quantities read
// their vectorizable parameter from the standard input and write to the standard
// output, so both streams are diverted for the duration of the call.
std::string query
(
    std::string const & quantity,
    std::map<std::string,std::string> const & params,
    std::string const & input = ""
)
{
    std::ostringstream out;
    std::istringstream in (input);

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

// Split the output table into its data rows, dropping the header comments. The
// direction triplets of 'tdcs' and 'ionamp' print as "( x y z )", so a row is
// returned as its whitespace-separated fields and the callers index from the end.
std::vector<std::vector<std::string>> rows (std::string const & table)
{
    std::vector<std::vector<std::string>> data;

    std::istringstream iss (table);
    std::string line;

    while (std::getline(iss, line))
    {
        std::istringstream fields (line);
        std::vector<std::string> row;
        std::string field;

        while (fields >> field)
            row.push_back(field);

        if (not row.empty() and row.front().front() != '#')
            data.push_back(row);
    }

    return data;
}

// value of the last-but-'back' field of the single data row of a table
double value (std::string const & table, std::size_t back = 0)
{
    std::vector<std::vector<std::string>> data = rows(table);

    if (data.size() != 1 or data.front().size() <= back)
    {
        std::cerr << "Expected one data row, got:\n" << table << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return std::stod(*(data.front().end() - 1 - back));
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

// --------------------------------------------------------------------------------- //

// composite Simpson rule over [a,b] with an even number of intervals
template <class Functor> double simpson (Functor f, double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++)
        sum += (i % 2 ? 4.0 : 2.0) * f(a + i * h);

    return sum * h / 3.0;
}

// the reference partial cross section, summed over the channels of one L
double reference_ics (int L, int S)
{
    auto integrand = [&](double beta) -> double
    {
        double sum = 0;

        for (Channel const & c : channels) if (c.L == L)
            sum += sqrabs(c.G(std::sin(beta)));

        return sum;
    };

    return 0.25 * (2*S + 1) * simpson(integrand, 0., 0.25 * special::constant::pi, 4096) / testKi;
}

// --------------------------------------------------------------------------------- //

int main (void)
{
    // the tables are compared against formulae written in Rydbergs, atomic units
    // of length and radians, so do not leave the choice to the defaults
    Eunits = eUnit_Ry;
    Lunits = lUnit_au;
    Aunits = aUnit_rad;

    hex_initialize("hex-ionization.db");
    hex_new();

    // Write the expansions the way hex-ecs does: generated on [0,kₘₐₓ] and stored as
    // a blob of their coefficients. Both spins carry the same radial amplitudes, so
    // that a wrong (2S+1)/4 weight cannot cancel against a different amplitude.
    std::ofstream sql ("batch-ionization.sql");
    sql << "BEGIN TRANSACTION;\n";

    for (Channel const & c : channels)
    {
        Chebyshev<double,Complex> cb;
        cb.generate([&](double k1) { return c.G(k1 / testKmax); }, 64, 0., testKmax);

        for (int S = 0; S <= 1; S++)
        {
            sql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                << testNi << "," << testLi << "," << testMi << ","
                << c.L << "," << S << "," << testEi << ","
                << c.l1 << "," << c.l2 << ","
                << cb.coeffs().toBlob().c_str() << ");\n";
        }
    }

    sql << "COMMIT;\n";
    sql.close();

    hex_import("batch-ionization.sql");
    hex_update();

    std::cout << "Ionization quantities at Ei = " << testEi << " Ry"
                 " (E_tot = " << testEtot << " Ry, k_max = " << testKmax << ")" << std::endl;

    // ------------------------------------------------------------------------- //
    // 1. the integral cross section itself

    double ics_total[2] = { 0., 0. };

    for (int S = 0; S <= 1; S++)
    for (int L = 0; L <= 2; L++)
    {
        std::map<std::string,std::string> p =
        {
            { "ni", str(testNi) }, { "li", str(testLi) },
            { "mi", str(testMi) }, { "nf", "0" }, { "lf", "0" },
            { "mf", "0" }, { "S", str(S) }, { "ell", str(L) },
            { "Ei", "-1" }   // negative energy asks for the stored values verbatim
        };

        double found = value(query("ics", p));
        ics_total[S] += found;

        check("ics       L = " + std::to_string(L) + ", S = " + std::to_string(S),
              found, reference_ics(L, S));
    }

    // ------------------------------------------------------------------------- //
    // 2. the single differential cross section with respect to the energy share

    auto reference_sdcsE = [&](double x, int S) -> double
    {
        double sum = 0;

        for (Channel const & c : channels)
            sum += sqrabs(c.G(std::sqrt(x)));

        double k1 = std::sqrt(testEtot * x);
        double k2 = std::sqrt(testEtot * (1 - x));

        return 0.25 * (2*S + 1) * sum * testEtot / (2 * k1 * k2 * testKi);
    };

    auto sdcsE = [&](double x, int S) -> double
    {
        std::map<std::string,std::string> p =
        {
            { "ni", str(testNi) }, { "li", str(testLi) },
            { "mi", str(testMi) }, { "S", str(S) },
            { "Ei", str(testEi) }, { "Eshare", str(x) }
        };

        return value(query("sdcsE", p));
    };

    for (int S = 0; S <= 1; S++)
    for (double x : { 0.1, 0.3, 0.5 })
    {
        check("sdcsE     x = " + std::to_string(x) + ", S = " + std::to_string(S),
              sdcsE(x, S), reference_sdcsE(x, S));
    }

    // ------------------------------------------------------------------------- //
    // 2b. the radial part of the amplitude on its own

    for (Channel const & c : channels)
    for (double x : { 0.2, 0.5 })
    {
        std::map<std::string,std::string> p =
        {
            { "ni", str(testNi) }, { "li", str(testLi) }, { "mi", str(testMi) },
            { "L", str(c.L) }, { "S", "0" },
            { "l1", str(c.l1) }, { "l2", str(c.l2) },
            { "Ei", str(testEi) }, { "Eshare", str(x) }
        };

        std::string out = query("ionf", p);
        double re = value(out, 1), im = value(out, 0);

        double k1 = std::sqrt(testEtot * x);
        double k2 = std::sqrt(testEtot * (1 - x));
        Complex expected = c.G(std::sqrt(x)) / (k1 * k2);

        check("ionf      (" + std::to_string(c.l1) + "," + std::to_string(c.l2) +
              ") x = " + std::to_string(x) + ", Re",
              re, expected.real());
        check("ionf      (" + std::to_string(c.l1) + "," + std::to_string(c.l2) +
              ") x = " + std::to_string(x) + ", Im",
              im, expected.imag());
    }

    // ------------------------------------------------------------------------- //
    // 3. integrating it over the energy share must give the integral cross section
    //
    // The substitution x = sin²β turns the integrable 1/√(x(1-x)) singularity of the
    // integrand into a constant. The first node is nudged off β = 0, where the query
    // would return an infinity that the vanishing Jacobian could not tame; the
    // integrand is even in β, so the displacement costs O(β²).

    for (int S = 0; S <= 1; S++)
    {
        auto integrand = [&](double beta) -> double
        {
            if (beta <= 0)
                beta = 1e-6;

            return sdcsE(std::sin(beta) * std::sin(beta), S) * std::sin(2 * beta);
        };

        check("∫ sdcsE dx           S = " + std::to_string(S),
              simpson(integrand, 0., 0.25 * special::constant::pi, 64), ics_total[S]);
    }

    // ------------------------------------------------------------------------- //
    // 4. and so must integrating sdcs12 over the cosine of the relative angle

    auto sdcs12 = [&](double theta, int S) -> double
    {
        std::map<std::string,std::string> p =
        {
            { "ni", str(testNi) }, { "li", str(testLi) },
            { "mi", str(testMi) }, { "S", str(S) },
            { "Ei", str(testEi) }, { "theta12", str(theta) }
        };

        return value(query("sdcs12", p));
    };

    for (int S = 0; S <= 1; S++)
    {
        auto integrand = [&](double c) -> double
        {
            return sdcs12(std::acos(std::max(-1.0, std::min(1.0, c))), S);
        };

        check("∫ sdcs12 d(cos θ₁₂)  S = " + std::to_string(S),
              simpson(integrand, -1., 1., 64), ics_total[S]);
    }

    // ------------------------------------------------------------------------- //
    // 5. Pinned sdcs12 values. Everything above averages the relative Coulomb
    //    phases away, because only λ = 0 survives the integration over the relative
    //    angle and its coefficient is diagonal in the angular channels. These two
    //    numbers were recorded from the fixed code; they move if the momenta that
    //    enter those phases, or the angular algebra around them, change.

    check("sdcs12    θ₁₂ = 1, S = 0", sdcs12(1.0, 0), 0.189948);
    check("sdcs12    θ₁₂ = 2, S = 1", sdcs12(2.0, 1), 0.897863);

    // ------------------------------------------------------------------------- //
    // 6. the triply differential cross section against the ionization amplitude

    for (int S = 0; S <= 1; S++)
    for (double share : { 0.25, 0.5 })
    {
        // ( θ₁ φ₁ E₁ ) ( θ₂ φ₂ E₂ ), the energies given as a ratio
        std::string dirs = "( 0.4 0.0 " + str(share) + " ) "
                           "( 2.1 1.3 " + str(1 - share) + " )\n";

        std::map<std::string,std::string> p =
        {
            { "ni", str(testNi) }, { "li", str(testLi) },
            { "mi", str(testMi) }, { "S", str(S) },
            { "Ei", str(testEi) }
        };

        // the amplitude is the last two fields, the cross section the third from last
        std::string amp = query("ionamp", p, dirs);
        double re = value(amp, 1), im = value(amp, 0);
        double tdcs = value(query("tdcs", p, dirs), 2);

        double k1 = std::sqrt(testEtot * share);
        double k2 = std::sqrt(testEtot * (1 - share));

        check("tdcs vs ionamp  share = " + std::to_string(share) + ", S = " + std::to_string(S),
              tdcs, 0.25 * (2*S + 1) * (re*re + im*im) * k1 * k2 / testKi);
    }

    // ------------------------------------------------------------------------- //

    if (failures > 0)
    {
        std::cout << failures << " check(s) failed." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "All ionization checks passed." << std::endl;
    return EXIT_SUCCESS;
}
