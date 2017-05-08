//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_VTKFILE_H
#define HEX_VTKFILE_H

// --------------------------------------------------------------------------------- //

#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"

// --------------------------------------------------------------------------------- //

class VTKRectGridFile
{
    public:
        
        void setGridX (rArray const & grid);
        void setGridY (rArray const & grid);
        void setGridZ (rArray const & grid);
        
        void appendScalarAttribute
        (
            std::string name,
            rArray const & array
        );
        
        void appendVector2DAttribute
        (
            std::string name,
            rArray const & xcomp,
            rArray const & ycomp
        );
        
        void appendVector3DAttribute
        (
            std::string name,
            rArray const & xcomp,
            rArray const & ycomp,
            rArray const & zcomp
        );
        
        void write (std::string filename) const;
    
    private:
        
        rArray xgrid_, ygrid_, zgrid_;
        std::vector<std::pair<std::string,std::vector<rArray>>> fields_;
};

// --------------------------------------------------------------------------------- //

#endif
