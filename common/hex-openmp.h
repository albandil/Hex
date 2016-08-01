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

#ifndef HEX_OPENMP
#define HEX_OPENMP

#include "hex-arrays.h"

// OpenMP macros
#ifdef _OPENMP
    #include <omp.h>
    
    #define OMP_CREATE_LOCKS(n) \
        std::vector<omp_lock_t> omplck(n); \
        for (omp_lock_t & L : omplck) \
            omp_init_lock(&L);

    #define OMP_LOCK_LOCK(i) \
        omp_set_lock(&omplck[i]);

    #define OMP_UNLOCK_LOCK(i) \
        omp_unset_lock(&omplck[i]);

    #define OMP_DELETE_LOCKS() \
        for (omp_lock_t & L : omplck) \
            omp_destroy_lock(&L);
#else
    #define OMP_CREATE_LOCKS(n)
    #define OMP_LOCK_LOCK(i)
    #define OMP_UNLOCK_LOCK(i)
    #define OMP_DELETE_LOCKS()
#endif

#endif
