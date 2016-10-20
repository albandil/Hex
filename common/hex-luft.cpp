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

#include "hex-luft.h"

// --------------------------------------------------------------------------------- //

std::shared_ptr<std::vector<LUft<LU_int_t,Complex>*>> LU;

// --------------------------------------------------------------------------------- //

LUftData defaultLUftData =
{
    /* drop_tolerance */            1e-8
    /* out_of_core */               , false
    /* verbosity */                 , 0
    /* Fortran MPI communicator */  , 0
    /* superlu_dist_grid */         , nullptr
};

// --------------------------------------------------------------------------------- //

template<> LUft<LU_int_t,Complex> * LUft<LU_int_t,Complex>::Choose (std::string factorizer)
{
    if (LU.get() == nullptr)
        HexException("No LU factorization method available.");
    
    for (LUft<LU_int_t,Complex> *lu : *LU)
    {
        if (lu->name() == factorizer or factorizer == "any")
        {
            return lu->New();
        }
    }
    
    std::cout << "Error!" << std::endl;
    std::cout << "  The LU factorizer \"" << factorizer << "\" is not available." << std::endl;
    std::cout << "  The program may not be compiled with support for that factorizer." << std::endl;
    std::cout << "  The available factorizers are: ";
    
    for (LUft<LU_int_t,Complex> *lu : *LU)
        std::cout << "\"" << lu->name() << "\" ";
    
    std::cout << std::endl << std::endl;
    
    std::exit(1);
}
