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

/**
 * @brief Write 3D dataset to a VTK file.
 * 
 * Write 3D grid data to a VTK file as a point dataset.
 * This file can be easily diplayed by the
 * free program ParaView. The expected storage of elements is
 * @f$ f(x_i, y_j, z_k) @f$ = f[i + j*nx + k*nx*ny],
 * where "nx", "ny" and "nz" are the lengths of the grids.
 */
void writeVTK_points
(
    std::ofstream & out,
    const ArrayView<Complex> f,
    const ArrayView<Real> xgrid,
    const ArrayView<Real> ygrid,
    const ArrayView<Real> zgrid,
    bool with_header = true
);

/**
 * @brief Write 3D dataset to a VTK file.
 * 
 * Write 3D grid data to a VTK file as a cell dataset.
 * This file can be easily diplayed by the
 * free program ParaView. The expected storage of elements is
 * @f$ f(\tilde{x}_i, \tilde{y}_j, \tilde{z}_k) @f$
 * = f[i + j*nx + k*nx*ny], where "nx", "ny" and "nz" are the
 * lengths of the grids and the tilded coordinates stand for the
 * centre of the [i,j,k]-th cell, not the [i,j,k]-th point (!)
 */
void writeVTK_cells
(
    std::ofstream & out,
    const ArrayView<Complex> f,
    const ArrayView<Real> xgrid,
    const ArrayView<Real> ygrid,
    const ArrayView<Real> zgrid
);

// --------------------------------------------------------------------------------- //

/**
 * @brief Another interface to writing VTK files.
 * 
 * This class can be used for the same purpose as the routines @ref writeVTK_cells
 * and @ref writeVTK_points.
 */
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
        
        void writePoints (std::string filename) const;
        void writeCells (std::string filename) const;
    
    private:
        
        rArray xgrid_, ygrid_, zgrid_;
        std::vector<std::pair<std::string,std::vector<rArray>>> fields_;
};

// --------------------------------------------------------------------------------- //

#endif
