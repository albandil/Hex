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

#ifndef HEX_DB_INTERFACES
#define HEX_DB_INTERFACES

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

/**
   @brief Initialize the environment.
   
   The function will open existing database with path specified in "dbname".
   If the database does not exist, a new empty database will be created.
   Also, all registered variables will be given chance to initialize by
   the function
   @code
      Variable::initialize (sqlitepp::session & db)
   @endcode
   This can be used to register special functions needed later in the SQL
   statements. For example the variable IntegralCrossSection registeres
   two functions: the simple square root and a more complicated numerical
   integrator that uses ClenshawCurtis quadrature.
   
   The function "initialize" also disables journalling to speed up insertion
   of data.
   
   @param dbname Zero-terminated character string.
*/
void hex_initialize_ (const char* dbname);
void hex_initialize (const char* dbname);

/**
   @brief Create new database.
   
   This will create the necessary table structure in a new empty database.
   It is necessary to call @ref initialize first, so that the databse file
   is created and opened. The tables are created by the call to
   @code
       Variable::SQL_CreateTable();
   @endcode
   for every registered variable. The return value is a vector of SQL statement
   which are executed.
*/
void hex_new_ ();
void hex_new ();

/**
   @brief Import SQL batch file.
   
   The effect of this function should be equivalent to direct use of sqlite3
   program:
   @code
       sqlite3 hex.db < batchfile.sql
   @endcode
   
   The function uses one line at a time and splits the lines into statements
   on semicolons. For this reason it does not allow multi-line statements.
   The statements produced by computational modules have always one statement
   per line.
   
   @param sqlname Zero-terminated character string.
*/
void hex_import_ (const char* sqlname);
void hex_import (const char* sqlname);

/**
   @brief Update integral and complete cross section.
   
   This function is used after import of T-matrices. All variables will be
   given chance to update themselves using the function
   @code
       Variable::SQL_Update();
   @endcode
   which returns SQL statement that should do the update. If a variable
   uses other variables besides T-matrices for its update, it must be inserted
   to VariableList AFTER those variables, because the order of update is
   set by the order in the list.
*/
void hex_update_ ();
void hex_update ();

/**
   @brief Optimize the SQLite database file.
   
   Reduces the occupied space using the VACUUM statement.
*/
void hex_optimize_ ();
void hex_optimize ();

/**
   @brief Dump contents of the T-matrix table.
   
   The output can be used to construct an equivalent table.
   The corresponding code is
   @code
   sqlite> .mode insert
   sqlite> select * from tmat
   @endcode
*/
void hex_dump_ (const char* dumpname);
void hex_dump (const char* dumpname);

/**
   @brief Scattering anplitude (Fortran).
   
   Fortran prototype equivalent to
   @code{.f90}
   subroutine scattering_amplitude (ni,li,mi,nf,lf,mf,S,E,N,angles,result)
     integer, intent(in) :: ni,li,mi
     integer, intent(in) :: nf,lf,mf
     integer, intent(int) :: S,N
     double precision, intent(in)     :: E
     double precision, dimension(N)   :: angles
     double precision, dimension(2*N) :: result
   @endcode
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param S Total spin (0 = singlet, 1 = triplet).
   @param nEnergies Number of supplied impact energies.
   @param energies Real array of length @c nEnergies containing impact energies in Rydbergs.
   @param nAngles Number of supplied angles.
   @param angles Real array of length @c nAngles containing scattering angles.
   @param result Complex array of length @c nEnergies*nAngles (or real array of twice that length) to contain the amplitudes.
   @param extra Extrapolated result (w.r.t. partial waves); optional - can be set to zero.
*/
void hex_scattering_amplitude_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S,
    int * nEnergies, double * energies,
    int * nAngles, double * angles,
    double * result, double * extra
);

/**
   @brief Scattering anplitude (C).
   
   C prototype.
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param S Total spin (0 = singlet, 1 = triplet).
   @param nEnergies Number of supplied impact energies.
   @param energies Real array of length @c nEnergies containing impact energies in Rydbergs.
   @param nAngles Number of supplied angles.
   @param angles Real array of length @c nAngles containing scattering angles.
   @param result Complex array of length @c nEnergies*nAngles (or real array of twice that length) to contain the amplitudes.
   @param extra Extrapolated result (w.r.t. partial waves); optional - can be set to zero.
*/
void hex_scattering_amplitude
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S,
    int nEnergies, double * energies,
    int nAngles, double * angles,
    double * result, double * extra
);

/**
   @brief Scattering anplitude for non-aligned impact direction (Fortran).
   
   This function will evaluate the scattering amplitude when the projectile is coming in a
   direction different from the quantization axis, i.e. when
   @f[
       \mathbf{k}_i \neq (0, 0, k_i) \,.
   @f]
   In that case the of the scattering amplitude can be computed using Wigner d-function
   and all possible amplitudes from scattering between various magnetic levels. The reason
   for this is that such situation is equivalent to a change of the quantization axis by the
   same amount. To be precise, it is
   @f[
       T_{n_f l_f m_f \leftarrow n_i l_i m_i} = \sum_{m_i' m_f'}
       D_{m_i' m_i}^{l_i} D_{m_f' m_f}^{l_f \ast} T_{n_f l_f m_f' \leftarrow n_i l_i m_i'} \,.
   @f]
   
   Fortran prototype equivalent to
   @code{.f90}
   subroutine scattering_amplitude (ni,li,mi,nf,lf,mf,S,E,N,angles,result)
     integer, intent(in) :: ni,li,mi
     integer, intent(in) :: nf,lf,mf
     integer, intent(int) :: S,N
     double precision, intent(in)     :: E
     double precision, dimension(N)   :: angles
     double precision, dimension(2*N) :: result
   @endcode
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param S Total spin (0 = singlet, 1 = triplet).
   @param E Impact energy in Rydbergs.
   @param N Sample count.
   @param alpha Impact angle (first of Euler angles).
   @param beta Impact angle (second of Euler angles).
   @param gamma Impact angle (third of Euler angles).
   @param angles Real array of length N containing scattering angles.
   @param result Complex array of length N (or real array of length 2N) to contain the amplitudes.
*/
void hex_scattering_amplitude_dir_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S, double * E, int * N,
    double * alpha, double * beta, double * gamma,
    double * angles, double * result
);

/**
   @brief Scattering anplitude for non-aligned impact direction (C).
   
   C prototype.
   
   See @ref hex_scattering_amplitude_dir_ for theory.
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param S Total spin (0 = singlet, 1 = triplet).
   @param E Impact energy in Rydbergs.
   @param N Sample count.
   @param alpha Impact angle (first of Euler angles).
   @param beta Impact angle (second of Euler angles).
   @param gamma Impact angle (third of Euler angles).
   @param angles Real array of length N containing scattering angles.
   @param result Complex array of length N (or real array of length 2N) to contain the amplitudes.
*/
void hex_scattering_amplitude_dir
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S, double E, int N,
    double alpha, double beta, double gamma,
    double * angles, double * result
);

/**
   @brief Differential cross section (Fortran).
   
   Fortran prototype equivalent to
   @code{.f90}
   subroutine hex_differential_cross_section (ni,li,mi,nf,lf,mf,S,E,N,angles,dcs)
     integer, intent(in)  :: ni,li,mi
     integer, intent(in)  :: nf,lf,mf
     integer, intent(in)  :: S,N
     double precision, intent(in)     :: E
     double precision, dimension(N)   :: angles
     double precision, dimension(2*N) :: dcs
   @endcode
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param S Total spin (0 = singlet, 1 = triplet).
   @param E Impact energy in Rydbergs.
   @param N Sample count.
   @param angles Real array of length N containing scattering angles.
   @param dcs Real array of length N to contain the cross sections.
*/
void hex_differential_cross_section_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * S,
    int * nEnergies, double * energies,
    int * nAngles, double * angles,
    double * dcs, double * extra
);

/**
   @brief Differential cross section (C).
   
   C prototype.
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param S Total spin (0 = singlet, 1 = triplet).
   @param E Impact energy in Rydbergs.
   @param N Sample count.
   @param angles Real array of length N containing scattering angles.
   @param dcs Real array of length N to contain the cross sections.
*/
void hex_differential_cross_section
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S,
    int nEnergies, double * energies,
    int nAngles, double * angles,
    double * dcs, double * extra
);

/**
   @brief Complete cross section (Fortran).
   
   Fortran prototype equivalent to
   @code{.f90}
   subroutine complete_cross_section (ni,li,mi,nf,lf,mf,N,energies,ccs,Nall)
     integer, intent(in)  :: ni,li,mi
     integer, intent(in)  :: nf,lf,mf
     integer, intent(in)  :: N
     double precision, dimension(N)   :: energies
     double precision, dimension(2*N) :: ccs
     integer, intent(out) :: Nall
   @endcode
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param N Sample count.
   @param energies Real array of impact energies in Rydbergs.
   @param ccs Real array of length N to contain the cross sections.
   @param Nall Number of energies available in database; will be set if @c ccs is set to null pointer.
   @param pws List of requested partial waves (or null pointer for all).
   @param npws Size of @c pws.
*/
void hex_complete_cross_section_
(
    int * ni, int * li, int * mi,
    int * nf, int * lf, int * mf,
    int * N, double * energies,
    double * ccs, int * Nall,
    int * npws, int * pws
);

/**
   @brief Complete cross section (C).
   
   C prototype.
   
   @param ni Initial atomic principal quantum number.
   @param li Initial atomic orbital quantum number.
   @param mi Initial atomic magnetic quantum number.
   @param nf Final atomic principal quantum number.
   @param lf Final atomic orbital quantum number.
   @param mf Final atomic magnetic quantum number.
   @param N Sample count.
   @param energies Real array of impact energies in Rydbergs.
   @param ccs Real array of length N to contain the cross sections.
   @param Nall Number of energies available in database; will be set if @c ccs is set to null pointer.
   @param pws List of requested partial waves (or null pointer for all).
   @param npws Size of @c pws.
*/
void hex_complete_cross_section
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int N, double * energies,
    double * ccs, int * Nall,
    int npws, int * pws
);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif
