/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
        std::cout << "\nUsage:\n\t./special f <lambda> <l1> <l2> <l1p> <l2p> <L>\n\n";
        return EXIT_FAILURE;
    }
    
    int lam = std::atoi(argv[2]);
    int l1  = std::atoi(argv[3]);
    int l2  = std::atoi(argv[4]);
    int l1p = std::atoi(argv[5]);
    int l2p = std::atoi(argv[6]);
    int L   = std::atoi(argv[7]);
    
    std::cout << "f[" << lam << "](" << l1 << "," << l2 << "," << l1p << ","
              << l2p << ") = " << computef(lam,l1,l2,l1p,l2p,L) << std::endl;
    
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
              << L << "," << M << "] = " << ClebschGordan(l1,m1,l2,m2,L,M) << std::endl;
    
    return EXIT_SUCCESS;
}

int run_Gaunt (int argc, char* argv[])
{
    if (argc < 8)
    {
        std::cout << "\nUsage:\n\t./special Gaunt <l1> <m1> <l2> <m2> <L> <M>\n\n";
        return EXIT_FAILURE;
    }
    
    int l1 = std::atoi(argv[2]);
    int m1 = std::atoi(argv[3]);
    int l2 = std::atoi(argv[4]);
    int m2 = std::atoi(argv[5]);
    int L  = std::atoi(argv[6]);
    int M  = std::atoi(argv[7]);
    
    std::cout << "G[" << l1 << "," << m1 << ";" << l2 << "," << m2 << ";"
              << L << "," << M << "] = " << Gaunt(l1,m1,l2,m2,L,M) << std::endl;
    
    return EXIT_SUCCESS;
}

int main (int argc, char* argv[])
{
    std::string name = (argc > 1 ? argv[1] : "help");
    bool needHelp = (std::find(help.begin(), help.end(), name) != help.end());
    
    if (needHelp)
    {
        std::cout << "\nEvaluates a special function as implemented in Hex.\n";
        std::cout << "\nUsage:\n\t./special <name> [options]\n\n";
        std::cout << "Available functions:\n";
    }
    
    Use(f);             // reduced matrix element factor
    Use(ClebschGordan); // Clebsch-Gordan coefficient
    Use(Gaunt);         // Gaunt coefficient
    
    if (needHelp)
    {
        std::cout << std::endl;
    }
    else
    {
        throw exception ("The function \"%s\" is not available.", argv[1]);
    }
    
    return EXIT_SUCCESS;
}
