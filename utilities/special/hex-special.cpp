//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include <cstdlib>
#include <iostream>

#include "hex-misc.h"
#include "hex-special.h"

std::vector<std::string> help = {
    std::string("-h"),
    std::string("--help"),
    std::string("-help"),
    std::string("help")
};

#define Stringify(x) #x
#define ToString(x) Stringify(x)
#define Use(fun) \
{ \
    if (name == ToString(fun)) \
        return run_##fun(argc,argv); \
    else if (needHelp) \
        std::cout << "\t* " << ToString(fun) << std::endl; \
}

int run_f (int argc, char* argv[])
{
    if (argc < 8)
    {
        std::cout << "\nUsage:\n\t./hex-special f <lambda> <l1> <l2> <l1p> <l2p> <L>\n\n";
        return EXIT_FAILURE;
    }

    int lam = std::atoi(argv[2]);
    int l1  = std::atoi(argv[3]);
    int l2  = std::atoi(argv[4]);
    int l1p = std::atoi(argv[5]);
    int l2p = std::atoi(argv[6]);
    int L   = std::atoi(argv[7]);

    std::cout << "f[" << lam << "](" << l1 << "," << l2 << "," << l1p << ","
              << l2p << ") = " << special::computef(lam,l1,l2,l1p,l2p,L) << std::endl;

    return EXIT_SUCCESS;
}

int run_ClebschGordan (int argc, char* argv[])
{
    if (argc < 8)
    {
        std::cout << "\nUsage:\n\t./special ClebschGordan <l1> <m1> <l2> <m2> <L> <M>\n\n";
        return EXIT_FAILURE;
    }

    int l1 = std::atoi(argv[2]);
    int m1 = std::atoi(argv[3]);
    int l2 = std::atoi(argv[4]);
    int m2 = std::atoi(argv[5]);
    int L  = std::atoi(argv[6]);
    int M  = std::atoi(argv[7]);

    std::cout << "C[" << l1 << "," << m1 << ";" << l2 << "," << m2 << ";"
              << L << "," << M << "] = " << special::ClebschGordan(l1,m1,l2,m2,L,M) << std::endl;

    return EXIT_SUCCESS;
}

int run_Gaunt (int argc, char* argv[])
{
    if (argc < 8)
    {
        std::cout << "\nUsage:\n\t./hex-special Gaunt <l1> <m1> <l2> <m2> <L> <M>\n\n";
        return EXIT_FAILURE;
    }

    int l1 = std::atoi(argv[2]);
    int m1 = std::atoi(argv[3]);
    int l2 = std::atoi(argv[4]);
    int m2 = std::atoi(argv[5]);
    int L  = std::atoi(argv[6]);
    int M  = std::atoi(argv[7]);

    std::cout << "G[" << l1 << "," << m1 << ";" << l2 << "," << m2 << ";"
              << L << "," << M << "] = " << special::Gaunt(l1,m1,l2,m2,L,M) << std::endl;

    return EXIT_SUCCESS;
}

int run_Y (int argc, char* argv[])
{
    if (argc < 6)
    {
        std::cout << "\nUsage:\n\t./hex-special Y [--degrees] <l> <m> <theta> <phi>\n\n";
        return EXIT_FAILURE;
    }

    int i = 2;
    double scale = 1;
    if (std::string(argv[i]) == std::string("--degrees"))
    {
        scale = special::constant::pi / 180.;
        i++;
    }

    if ((scale == 1. and argc != 6) or (scale != 1. and argc != 7))
    {
        std::cout << "\nUsage:\n\t./hex-special Y [--degrees] <l> <m> <theta> <phi>\n\n";
        return EXIT_FAILURE;
    }

    int l = std::atoi(argv[i++]);
    int m = std::atoi(argv[i++]);
    double theta = std::atof(argv[i++]);
    double phi   = std::atof(argv[i++]);

    std::cout << "Y[" << l << "," << m << "](" << theta << (scale == 1. ? "," : "°,")
              << phi << (scale == 1. ? ") = " : "°) = ")
              << special::sphY(l, m, theta*scale, phi*scale) << std::endl;

    return EXIT_SUCCESS;
}

int run_Coulomb (int argc, char* argv[])
{
    if (argc < 6)
    {
        std::cout << std::endl << "Usage:" << std::endl;
        std::cout << "\t./hex-special Coulomb --eta-rho <L> <eta> <rho>" << std::endl;
        std::cout << "\t./hex-special Coulomb --k-r <L> <k> <r>" << std::endl;
        std::cout << std::endl;
        return EXIT_FAILURE;
    }

    if (std::string(argv[2]) == std::string("--eta-rho"))
    {    
        int L = std::atoi(argv[3]);
        double eta = std::atof(argv[4]);
        double rho = std::atof(argv[5]);

        double F, expF;
        gsl_sf_coulomb_wave_F_array(L, 0, eta, rho, &F, &expF);

        std::cout << "F[" << L << "](" << eta << "," << rho << ") = " << F << std::endl;
    }
    else if (std::string(argv[2]) == std::string("--k-r"))
    {
        int L = std::atoi(argv[3]);
        double k = std::atof(argv[4]);
        double r = std::atof(argv[5]);

        double F, expF;
        gsl_sf_coulomb_wave_F_array(L, 0, -1/k, k*r, &F, &expF);

        std::cout << "F[" << L << "](" << k << "," << r << ") = " << F << std::endl;
    }
    else
    {
        std::cout << std::endl << "Usage:" << std::endl;
        std::cout << "\t./hex-special Coulomb --eta-rho <L> <eta> <rho>" << std::endl;
        std::cout << "\t./hex-special Coulomb --k-r <L> <k> <r>" << std::endl;
        std::cout << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int run_Wigner9j (int argc, char* argv[])
{
    if (argc < 11)
    {
        std::cout << std::endl << "Usage:" << std::endl;
        std::cout << "\t./hex-special Wigner9j <2*ja> <2*jb> <2*jc> <2*jd> <2*je> <2*jf> <2*jg> <2*jh> <2*ji>" << std::endl << std::endl;
        return EXIT_FAILURE;
    }

    int two_ja = std::atoi(argv[2]), two_jb = std::atoi(argv[3]), two_jc = std::atoi(argv[4]);
    int two_jd = std::atoi(argv[5]), two_je = std::atoi(argv[6]), two_jf = std::atoi(argv[7]);
    int two_jg = std::atoi(argv[8]), two_jh = std::atoi(argv[9]), two_ji = std::atoi(argv[10]);

    double w9j = gsl_sf_coupling_9j
    (
        two_ja, two_jb, two_jc,
        two_jd, two_je, two_jf,
        two_jg, two_jh, two_ji
    );

    std::cout << format
    (
        "Wigner9j(%g,%g,%g; %g,%g,%g; %g,%g,%g) = %g",
        0.5 * two_ja, 0.5 * two_jb, 0.5 * two_jc,
        0.5 * two_jd, 0.5 * two_je, 0.5 * two_jf,
        0.5 * two_jg, 0.5 * two_jh, 0.5 * two_ji,
        w9j
    ) << std::endl;

    return EXIT_SUCCESS;
}

int main (int argc, char* argv[])
{
    std::string name = (argc > 1 ? argv[1] : "help");
    bool needHelp = (std::find(help.begin(), help.end(), name) != help.end());

    if (needHelp)
    {
        std::cout << std::endl << "Evaluates a special function as implemented in Hex." << std::endl;
        std::cout << std::endl << "Usage:" << std::endl << "\t./hex-special <name> [options]" << std::endl;
        std::cout << std::endl << "Available functions:" << std::endl;
    }

    Use(f);             // reduced matrix element factor
    Use(ClebschGordan); // Clebsch-Gordan coefficient
    Use(Gaunt);         // Gaunt coefficient
    Use(Y);             // spherical harmonic
    Use(Coulomb);       // Coulomb function
    Use(Wigner9j);      // Wigner 9j-coefficient

    if (needHelp)
    {
        std::cout << std::endl;
        std::cout << "By calling \"./hex-special <name>\" (i.e. without further arguments) "
                      "you will receive detailed information on the chosen special function.";
        std::cout << std::endl << std::endl;
    }
    else
    {
        HexException("The function \"%s\" is not available.", argv[1]);
    }

    return EXIT_SUCCESS;
}
