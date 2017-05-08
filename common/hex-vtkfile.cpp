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

#include <fstream>
#include <iostream>

// --------------------------------------------------------------------------------- //

#include "hex-vtkfile.h"

// --------------------------------------------------------------------------------- //

void VTKRectGridFile::setGridX (rArray const & grid)
{
    xgrid_ = grid;
}

void VTKRectGridFile::setGridY (rArray const & grid)
{
    ygrid_ = grid;
}

void VTKRectGridFile::setGridZ (rArray const & grid)
{
    zgrid_ = grid;
}

void VTKRectGridFile::appendScalarAttribute (std::string name, rArray const & array)
{
    fields_.push_back(std::make_pair(name, std::vector<rArray>{ array }));
}

void VTKRectGridFile::appendVector2DAttribute (std::string name, rArray const & xcomp, rArray const & ycomp)
{
    fields_.push_back(std::make_pair(name, std::vector<rArray>{ xcomp, ycomp }));
}

void VTKRectGridFile::appendVector3DAttribute (std::string name, rArray const & xcomp, rArray const & ycomp, rArray const & zcomp)
{
    fields_.push_back(std::make_pair(name, std::vector<rArray>{ xcomp, ycomp, zcomp }));
}

void VTKRectGridFile::write (std::string filename) const
{
    std::ofstream out (filename);
    
    // write VTK header
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Wave function" << std::endl;
    out << "ASCII" << std::endl;
    
    // get sizes
    std::size_t nx = xgrid_.size();
    std::size_t ny = ygrid_.size();
    std::size_t nz = zgrid_.size();
    std::size_t N = nx * ny * nz;
    
    // write grid data
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    out << "X_COORDINATES " << nx << " float" << std::endl;
    out << to_string(xgrid_) << std::endl;
    out << "Y_COORDINATES " << ny << " float" << std::endl;
    out << to_string(ygrid_) << std::endl;
    out << "Z_COORDINATES " << nz << " float" << std::endl;
    out << to_string(zgrid_) << std::endl;
    out << "POINT_DATA " << N << std::endl;
    out << "FIELD wavefunction " << fields_.size() << std::endl;
    
    // write field data
    for (std::pair<std::string,std::vector<rArray>> const & field : fields_)
    {
        std::string const & fieldname = field.first;
        std::size_t ncomp = field.second.size();
        
        out << fieldname << " " << ncomp << " " << N << " float" << std::endl;
        for (unsigned i = 0; i < nx; i++)
        {
            for (unsigned j = 0; j < ny; j++)
            for (unsigned k = 0; k < nz; k++)
            {
                for (unsigned icomp = 0; icomp < ncomp; icomp++)
                    out << field.second[icomp][(i * ny + j) * nz + k] << " ";
            }
            out << std::endl;
        }
    }
}
