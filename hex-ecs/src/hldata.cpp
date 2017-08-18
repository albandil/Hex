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

#include "hldata.h"

void HlData::hdflink (const char* file)
{
    filename = file;
}

bool HlData::hdfcheck (const char* file) const
{
    // open HDF file for reading
    HDFFile hdf ((file == nullptr ? filename.c_str() : file), HDFFile::readonly);
    return hdf.valid();
}

bool HlData::hdfload (const char* file)
{
    // open HDF file for reading
    HDFFile hdf ((file == nullptr ? filename.c_str() : file), HDFFile::readonly);
    if (not hdf.valid())
        return false;
    
    // read size
    unsigned size;
    if (not hdf.read("n", &size, 1))
        return false;
    
    // allocate memory (if needed)
    if (Cl.rows() != (int)size or Cl.cols() != (int)size)
    {
        Cl.drop();
        Cl = std::move(ColMatrix<Complex>(size,size));
    }
    if (invCl_invsqrtS.rows() != (int)size or invCl_invsqrtS.cols() != (int)size)
    {
        invCl_invsqrtS.drop();
        invCl_invsqrtS = std::move(RowMatrix<Complex>(size,size));
    }
    if (invsqrtS_Cl.rows() != (int)size or invsqrtS_Cl.cols() != (int)size)
    {
        invsqrtS_Cl.drop();
        invsqrtS_Cl = std::move(RowMatrix<Complex>(size,size));
    }
    
    // read matrices
    if (not hdf.read("Cl", Cl.data().data(), size * size))
        return false;
    if (not hdf.read("invCl_invsqrtS", invCl_invsqrtS.data().data(), size * size))
        return false;
    if (not hdf.read("invsqrtS_Cl", invsqrtS_Cl.data().data(), size * size))
        return false;
    
    // read eigenvalues
    Dl.resize(size);
    if (not hdf.read("Dl", Dl.data(), size))
    {
        std::cout << "Failed to read Dl from " << hdf.name() << ": " << hdf.error() << std::endl;
        return false;
    }
    
    return true;
}

bool HlData::hdfsave (const char* file) const
{
    // open HDF file for writing
    HDFFile hdf ((file == nullptr ? filename.c_str() : file), HDFFile::overwrite);
    if (not hdf.valid())
        return false;
    
    // write size
    std::size_t size = Dl.size();
    if (not hdf.write("n", &size, 1))
        return false;
    
    // write matrices
    if (not hdf.write("Cl", Cl.data().data(), size * size))
        return false;
    if (not hdf.write("invCl_invsqrtS", invCl_invsqrtS.data().data(), size * size))
        return false;
    if (not hdf.write("invsqrtS_Cl", invsqrtS_Cl.data().data(), size * size))
        return false;
    
    // write eigenvalues
    if (not hdf.write("Dl", Dl.data(), size))
        return false;
    
    return true;
}

void HlData::drop ()
{
    invCl_invsqrtS.drop();
    invsqrtS_Cl.drop();
    Dl.clear();
}

cArray HlData::readPseudoState (unsigned l, unsigned ichan) const
{
    HDFFile datafile (filename, HDFFile::readonly);
    if (not datafile.valid())
        HexException("File %s with the one-electron eigenstates was not found. Run the solver again to regenerate it.", filename.c_str());
    
    // get rank of the hamiltonian
    std::size_t n;
    if (not datafile.read("n", &n, 1))
        HexException("File %s does not contain the requested dataset \"n\".", filename.c_str());
    if (ichan >= n)
        HexException("File %s contains only %d (<= %d) channels.", filename.c_str(), n, ichan);
    
    // read the eigenenergies
    cArray energies (n);
    if (not datafile.read("Dl", &energies[0], n))
        HexException("File %s does not contain the requested dataset \"Dl\".", filename.c_str());
    
    // sort energies
    iArray indices (n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j){ return energies[i].real() < energies[j].real(); });
    
    // read the right eigenstate
    cArray data (n);
    if (not datafile.read("Cl", &data[0], n, indices[ichan] * n))
        HexException("Failed to read the pseudostate l = %d, ichan = %d from the dataset \"Cl\" in file %s.", l, ichan, filename.c_str());
    
    // Adjust the overall sign of the eigenvector so that the result is compatible with the
    // sign convention of GSL's function gsl_sf_hydrogenicR (used in previous versions of hex-ecs).
    // That is, the radial function should increase from origin to positive values, then turn back
    // and (potentially) dive through zero.
    
    if (data.front().real() < 0.0_r)
    {
        for (Complex & z : data)
            z = -z;
    }
    
    return data;
}
