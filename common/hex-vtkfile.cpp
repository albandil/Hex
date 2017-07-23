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

void writeVTK_points
(
    std::ofstream & out,
    const cArrayView ev,
    const rArrayView xgrid,
    const rArrayView ygrid,
    const rArrayView zgrid,
    bool with_header
)
{
    // array lengths
    unsigned nx = xgrid.size();
    unsigned ny = ygrid.size();
    unsigned nz = zgrid.size();
    unsigned N = nx * ny * nz;
    
    // write VTK header
    if (with_header)
    {
        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "Wave function" << std::endl;
        out << "ASCII" << std::endl;
    }
    
    // write grid data
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    out << "X_COORDINATES " << nx << " float" << std::endl;
    out << to_string(xgrid) << std::endl;
    out << "Y_COORDINATES " << ny << " float" << std::endl;
    out << to_string(ygrid) << std::endl;
    out << "Z_COORDINATES " << nz << " float" << std::endl;
    out << to_string(zgrid) << std::endl;
    out << "POINT_DATA " << N << std::endl;
    out << "FIELD wavefunction " << 2 * ev.size() / N << std::endl;
    
    for (std::size_t iblock = 0; (iblock + 1) * N <= ev.size(); iblock++)
    {
        // save real part
        out << "block_" << iblock << "_re 1 " << N << " float" << std::endl;
        for (unsigned i = 0; i < nx; i++)
        {
            for (unsigned j = 0; j < ny; j++)
            for (unsigned k = 0; k < nz; k++)
                out << ev[((iblock * nx + i) * ny + j) * nz + k].real() << " ";
            out << std::endl;
        }
        
        // save imaginary part
        out << "block_" << iblock << "_im 1 " << N << " float" << std::endl;
        for (unsigned i = 0; i < nx; i++)
        {
            for (unsigned j = 0; j < ny; j++)
            for (unsigned k = 0; k < nz; k++)
                out << ev[((iblock * nx + i) * ny + j) * nz + k].imag() << " ";
            out << std::endl;
        }
    }
}

void writeVTK_cells
(
    std::ofstream & out,
    const cArrayView ev,
    const rArrayView xgrid,
    const rArrayView ygrid,
    const rArrayView zgrid
)
{
    // point counts, cell counts
    unsigned px = xgrid.size(), cx = px - 1;
    unsigned py = ygrid.size(), cy = py - 1;
    unsigned pz = zgrid.size(), cz = pz - 1;
    unsigned N = cx * cy * cz;
    
    assert (ev.size() == N);
    
    // write VTK header
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Wave function" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << px << " " << py << " " << pz << std::endl;
    out << "X_COORDINATES " << px << " float" << std::endl;
    out << to_string(xgrid) << std::endl;
    out << "Y_COORDINATES " << py << " float" << std::endl;
    out << to_string(ygrid) << std::endl;
    out << "Z_COORDINATES " << pz << " float" << std::endl;
    out << to_string(zgrid) << std::endl;
    out << "CELL_DATA " << N << std::endl;
    out << "FIELD wavefunction 2" << std::endl;
    
    // save real part
    out << "Re 1 " << N << " float" << std::endl;
    for (unsigned i = 0; i < cx; i++)
    {
        for (unsigned j = 0; j < cy; j++)
        for (unsigned k = 0; k < cz; k++)
            out << ev[(i * cy + j) * cz + k].real() << " ";
        out << std::endl;
    }
    
    // save imaginary part
    out << "Im 1 " << N << " float" << std::endl;
    for (unsigned i = 0; i < cx; i++)
    {
        for (unsigned j = 0; j < cy; j++)
        for (unsigned k = 0; k < cz; k++)
            out << ev[(i * cy + j) * cz + k].imag() << " ";
        out << std::endl;
    }
}

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

void VTKRectGridFile::writePoints (std::string filename) const
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

void VTKRectGridFile::writeCells (std::string filename) const
{
    std::ofstream out (filename);
    
    // write VTK header
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Wave function" << std::endl;
    out << "ASCII" << std::endl;
    
    // get number of points
    std::size_t nx = xgrid_.size();
    std::size_t ny = ygrid_.size();
    std::size_t nz = zgrid_.size();
    
    if (nx == 0 or ny == 0 or nz == 0)
        HexException("All grids must be non-empty");
    
    // get number of cells
    std::size_t cx = std::max<std::size_t>(nx - 1, 1);
    std::size_t cy = std::max<std::size_t>(ny - 1, 1);
    std::size_t cz = std::max<std::size_t>(nz - 1, 1);
    std::size_t Nc = cx * cy * cz;
    
    // write grid data
    out << "DATASET RECTILINEAR_GRID" << std::endl;
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    out << "X_COORDINATES " << nx << " float" << std::endl;
    out << to_string(xgrid_) << std::endl;
    out << "Y_COORDINATES " << ny << " float" << std::endl;
    out << to_string(ygrid_) << std::endl;
    out << "Z_COORDINATES " << nz << " float" << std::endl;
    out << to_string(zgrid_) << std::endl;
    out << "CELL_DATA " << Nc << std::endl;
    out << "FIELD wavefunction " << fields_.size() << std::endl;
    
    // write field data
    for (std::pair<std::string,std::vector<rArray>> const & field : fields_)
    {
        std::string const & fieldname = field.first;
        std::size_t ncomp = field.second.size();
        
        // check that all components have the right dimension
        for (unsigned icomp = 0; icomp < ncomp; icomp++)
        {
            if (field.second[icomp].size() != Nc)
            {
                HexException
                (
                    "Size of field %s (component %d) is not compatible with given bases (%d != %d x %d x %d).",
                    fieldname.c_str(), icomp, field.second[icomp].size(), cx, cy, cz
                );
            }
        }
        
        // write data to file
        out << fieldname << " " << ncomp << " " << Nc << " float" << std::endl;
        for (unsigned k = 0; k < cz; k++)
        {
            for (unsigned j = 0; j < cy; j++)
            for (unsigned i = 0; i < cx; i++)
            {
                for (unsigned icomp = 0; icomp < ncomp; icomp++)
                    out << field.second[icomp][(i * cy + j) * cz + k] << " ";
            }
            out << std::endl;
        }
    }
}
