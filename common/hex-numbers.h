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

#ifndef HEX_NUMBERS
#define HEX_NUMBERS

#include <cassert>
#include <complex>
#include <cmath>

#undef I

#ifdef WITH_MPI
#include <mpi.h>
#endif

//
// Choose integer type for use in LU-factorization libraries.
//

#ifndef _LONGINT
    // use 32-bit integers (supported by both UMFPACK and SuperLU)
    typedef std::int32_t LU_int_t;
#else
    // use 64-bit integers (not supported by serial SuperLU)
    typedef std::int64_t LU_int_t;
#endif

//
// Complex data type alias.
//

// shorthand for std::complex<double>
typedef std::complex<double> Complex;

// complex number literal
inline Complex operator "" _z (long double x)
{
    return Complex(x,0.);
}

// imaginary number literal
inline Complex operator "" _i (long double x)
{
    return Complex(0.,x);
}

/// Squared modulus of a real number.
inline double sqrabs (double z)
{
    return z*z;
}

/// Squared modulus of a complex number.
inline double sqrabs (Complex z)
{
    return z.real() * z.real() + z.imag() * z.imag();
}

/// Complex ordering by real parts.
inline bool Complex_realpart_less (Complex const & a, Complex const & b)
{
    return a.real() < b.real();
}

/// Complex ordering by imaginary parts.
inline bool Complex_imagpart_less (Complex const & a, Complex const & b)
{
    return a.imag() < b.imag();
}

/// Finite check for complex number.
inline bool Complex_finite (Complex const & z)
{
    return std::isfinite(z.real()) and std::isfinite(z.imag());
}

//
// Data types.
//

/**
 * @brief HDF data types.
 * 
 * These data types are used for Hex data format. The identificators are always a multiple
 * of 0x100 (= 256) plus the bytesize of the data type.
 */
typedef enum
{
    DataInt32    = 0x100 + sizeof(int),
    DataInt64    = 0x200 + sizeof(std::int64_t),
    DataUInt32   = 0x300 + sizeof(unsigned),
    DataUInt64   = 0x400 + sizeof(std::uint64_t),
    DataFloat32  = 0x500 + sizeof(float),
    DataDouble64 = 0x600 + sizeof(double),
    DataDouble80 = 0x700 + sizeof(long double)
}
HDFDataType;

//
// Data type abstract traits.
//

/**
 * @brief Data-type information.
 * 
 * The information about components of a specific data type can be used in
 * type-generic template functions. This class (or rather its specializations
 * for individual types) offer the necessary information.
 * 
 * An example of use would be a formatted output of column data of an array
 * of the abstract data type 'T'. If we wanted to print every component as a
 * separate column, we would need to (a) know the total number of components
 * of the data type 'T' and (b) have the tool to access individual elements.
 * This class (or rather its specializations for different types) offers
 * both through the members 'ncmpt' and 'cmpt'.
 */
template <class T> class typeinfo {};

/// Data-type info class specialization for 'int'.
template<> class typeinfo<int>
{
    public:
        /// Component data type.
        typedef int cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, int x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataInt32; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_INT; }
#endif
};

/// Data-type info class specialization for 'int64'.
template<> class typeinfo<std::int64_t>
{
    public:
        /// Component data type.
        typedef std::int64_t cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, int x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataInt64; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_INT64_T; }
#endif
};

/// Data-type info class specialization for 'unsigned int'.
template<> class typeinfo<unsigned>
{
    public:
        /// Component data type.
        typedef unsigned cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, int x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataUInt32; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_UNSIGNED; }
#endif
};

/// Data-type info class specialization for 'unsigned int64'.
template<> class typeinfo<std::uint64_t>
{
    public:
        /// Component data type.
        typedef std::uint64_t cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, int x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataUInt64; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_UINT64_T; }
#endif
};

/// Data-type info class specialization for 'float'.
template<> class typeinfo<float>
{
    public:
        /// Component data type.
        typedef float cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, float x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataFloat32; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_FLOAT; }
#endif
};

/// Data-type info class specialization for 'double'.
template<> class typeinfo<double>
{
    public:
        /// Component data type.
        typedef double cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, double x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataDouble64; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_DOUBLE; }
#endif
};

/// Data-type info class specialization for 'long double'.
template<> class typeinfo<long double>
{
    public:
        /// Component data type.
        typedef long double cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 1;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, double x) { assert(i < ncmpt); return x; }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return DataDouble80; }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return MPI_LONG_DOUBLE; }
#endif
};

/// Data-type info class specialization for 'std::complex'.
template<> template<class T> class typeinfo<std::complex<T>>
{
    public:
        /// Component data type.
        typedef T cmpttype;
        
        /// Component count.
        static const std::size_t ncmpt = 2;
        
        /// Component getter.
        static cmpttype cmpt (std::size_t i, std::complex<T> x) { assert(i < ncmpt); return (i == 0 ? x.real() : x.imag()); }
        
        /// HDF data type identification.
        static HDFDataType hdfcmpttype () { return typeinfo<T>::hdfcmpttype(); }
        
#ifdef WITH_MPI
        /// MPI data type of a component.
        static MPI_Datatype mpicmpttype () { return typeinfo<T>::mpicmpttype(); }
#endif
};

#endif
