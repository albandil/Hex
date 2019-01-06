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

#include <iostream>
#include <cstdlib>
#include <string>

#include "hex-arrays.h"
#include "hex-misc.h"
#include "hex-compact.h"
#include "hex-chebyshev.h"
#include "hex-hdffile.h"

template<typename T> void load_and_write(const char* hdf, int samples)
{
    // load coefficients
    NumberArray<T> coefs;
    if (not coefs.hdfload(hdf))
    {
        std::cerr << "Can't read file \"" << hdf << "\"\n";
        std::exit(EXIT_FAILURE);
    }

    // create the expansion
    Chebyshev<double,T> expansion(coefs, -1, 1);

    // evaluate the expansion
    for (int i = 0; i <= samples; i++)
    {
        double x = (2.*i-samples) / samples;

        Complex y = expansion.clenshaw(x,coefs.size());
        std::cout << x << "\t" << y.real() << "\t" << y.imag() << "\n";
    }
}

int main (int argc, char *argv[])
{
    bool cpx = false;    // whether to zip complex expansion
    std::string hdf;
    int samples = -1;

    for (int iarg = 1; iarg < argc; iarg++)
    {
        if (strcmp(argv[iarg],"--complex") == 0)
        {
            cpx = true;
        }
        else if (hdf.empty())
        {
            hdf = std::string(argv[iarg]);
        }
        else
        {
            char* tail;
            samples = strtol(argv[iarg], &tail, 10);
            if (*tail != 0)
                samples = -1;
        }
    }

    if (samples < 0)
    {
        std::cerr << "\nUsage: ./chebeval [--complex] <HDFfile> <samples>\n\n";
        std::exit(EXIT_FAILURE);
    }

    if (cpx)
        load_and_write<Complex>(hdf.c_str(), samples);
    else
        load_and_write<double>(hdf.c_str(), samples);

    return 0;
}
