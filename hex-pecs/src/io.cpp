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

#include "hex-special.h"

#include "io.h"

CommandLine::CommandLine (int argc, char* argv[])
{
    forward_grid = true;
    forward_states = true;
    backward = true;
    inputfile = "pecs.inp";
}

InputFile::InputFile(const CommandLine& cmd)
{
    istates.push_back(HState({ 1, 0, 0 }));
    fstates.push_back(HState({ 1, 0, 0 }));
    Etot = 1;
    Z = 1;
    L = Pi = S = nL = 0;
    ecstheta = special::constant::pi_quart;
    rgrid = concatenate
    (
//         linspace(0.00, 0.20, 21),
//         linspace(0.26, 2.00, 30),
//         linspace(2.20, 20.0, 90),
//         linspace(20.4, 200., 450)
        linspace(0.,3.,4)
    );
    cgrid = concatenate
    (
//         linspace(0.1, 2.0, 20),
//         linspace(2.5, 20.0, 36)
        linspace(0.,3.,4)
    );
}
