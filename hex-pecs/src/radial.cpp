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

#include "radial.h"

RadialBasis::RadialBasis (InputFile const & inp)

{
    Complex rotation (std::cos(inp.ecstheta), std::sin(inp.ecstheta));
    
    rgrid = inp.rgrid;
    cgrid = inp.cgrid * rotation;
    
    std::cout << "Complex rotation: angle = " << inp.ecstheta << ", factor = " << rotation << std::endl;
    
    Npts = rgrid.size() + cgrid.size() - 1;
    grid.resize(Npts);
    
    for (std::size_t i = 0; i < Npts; i++)
    {
        if (i < rgrid.size())
            grid[i] = rgrid[i];
        else
            grid[i] = cgrid[i - rgrid.size() + 1] + rgrid.back();
    }
    
    std::cout << "Real grid points: " << std::endl;
    for (std::string line : inp.rgrid.lines(100))
        std::cout << '\t' << line << std::endl;
    std::cout << std::endl;
    
    std::cout << "Complex grid points: " << std::endl;
    for (std::string line : (inp.cgrid + inp.rgrid.back()).lines(100))
        std::cout << '\t' << line << std::endl;
    std::cout << std::endl;
    
    std::cout << "Full grid points: " << std::endl;
    for (std::string line : grid.lines(100))
        std::cout << '\t' << line << std::endl;
    std::cout << std::endl;
}
