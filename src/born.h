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

#ifndef HEX_BORN_SYMB
#define HEX_BORN_SYMB

#include <map>

#include <ginac/ginac.h>

#include "numbers.h"
#include "vec3d.h"

GiNaC::ex rl_Y
(
    int l, int m,
    GiNaC::symbol const & rp, GiNaC::symbol const & rm,
    GiNaC::realsymbol const & rz
);

GiNaC::ex psi_nlm_poly
(
    int n, int l, int m,
    GiNaC::possymbol const & r,
    GiNaC::symbol const & rp, GiNaC::symbol const & rm,
    GiNaC::realsymbol const & rz
);

std::map<std::tuple<int,int,int,int,int>,Complex> Wb_symb_in
(
    int n1, int l1, int m1,
    int n2, int l2, int m2
);

Complex eval_Wb
(
    std::map<std::tuple<int,int,int,int,int>,Complex> const & poly,
    double nu, geom::vec3d const & k
);

std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> W_symb_in
(
    int n, int l, int m
);

Complex eval_W
(
    std::map<std::tuple<int,int,int,int,int,int,int,int,int,int>,Complex> const & poly,
    double nu,
    geom::vec3d const & vk,
    geom::vec3d const & vq
);

#endif // HEX_BORN_SYMB
