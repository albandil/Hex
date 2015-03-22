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

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/hdffile.cpp"

int main (int argc, char *argv[])
{
    bool cpx = false;    // whether to write complex expansion
    std::string hdf;
    
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
            break;
    }
    
    if (argc < 2 or hdf.empty())
    {
        std::cout << "\nUsage:\n\t./hdf2txt [--complex] <HDFfile>\n\n";
        std::exit(EXIT_SUCCESS);
    }
    
    if (cpx)
    {
        // read HDF file
        cArray a;
        a.hdfload(hdf);
    
        // write raw data
        for (int i = 0; i < a.size(); i++)
            std::cout << a[i] << "\n";
    }
    else
    {
        // read HDF file
        rArray a;
        a.hdfload(hdf);
    
        // write raw data
        for (int i = 0; i < a.size(); i++)
            std::cout << a[i] << "\n";
    }
    
    return 0;
}
