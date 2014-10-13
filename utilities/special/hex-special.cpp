//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#include "../../src/special.h"
#include "../../src/special.cpp"

#include <cstdlib>
#include <iostream>

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

int main (int argc, char* argv[])
{
    std::string name = (argc > 1 ? argv[1] : "help");
    bool needHelp = (std::find(help.begin(), help.end(), name) != help.end());
    
    if (needHelp)
    {
        std::cout << "\nEvaluates a special function as implemented in Hex.\n";
        std::cout << "\nUsage:\n\t./hex-special <name> [options]\n\n";
        std::cout << "Available functions:\n";
    }
    
    Use(f);             // reduced matrix element factor
    Use(ClebschGordan); // Clebsch-Gordan coefficient
    Use(Gaunt);         // Gaunt coefficient
    Use(Y);             // spherical harmonic
    Use(Coulomb)        // Coulomb function
    
    if (needHelp)
    {
        std::cout << "By calling \"./hex-special <name>\" (i.e. without further arguments) "
                      "you will receive detailed information on the chosen special function."
                  << std::endl << std::endl;
    }
    else
    {
        throw exception ("The function \"%s\" is not available.", argv[1]);
    }
    
    return EXIT_SUCCESS;
}
