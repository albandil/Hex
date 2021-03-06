//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2018, Jakub Benda, Charles University in Prague                    //
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

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"
#include "hex-chebyshev.h"
#include "hex-clenshawcurtis.h"
#include "hex-hdffile.h"
#include "hex-hydrogen.h"
#include "hex-matrix.h"
#include "hex-misc.h"
#include "hex-special.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "amplitudes.h"
#include "bspline.h"
#include "hldata.h"
#include "radial.h"

// --------------------------------------------------------------------------------- //

// Extrapolate uniformly sampled function y(x) = a/x + b.
Complex inv_power_extrapolate (rArray X, cArrayView Y)
{
    // invert independent variable
    for (Real & x : X)
        x = 1.0 / x;

    // do the regression
    Real avg_x   = sum(X)   / X.size();
    Real avg_xx  = sum(X*X) / X.size();
    Complex avg_y  = sum(Y)   / Real(X.size());
    Complex avg_xy = sum(X*Y) / Real(X.size());
    Complex a = (avg_xy - avg_x * avg_y) / (avg_xx - avg_x * avg_x);
    Complex b = avg_y - a * avg_x;

    // return the limit
    return b;
}

// --------------------------------------------------------------------------------- //

Amplitudes::Amplitudes
(
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    InputFile const & inp,
    Parallel const & par,
    CommandLine const & cmd,
    std::vector<std::pair<int,int>> const & ang
) : bspline_inner_(bspline_inner), bspline_full_(bspline_full),
    rad_(bspline_inner, bspline_full, 0),
    inp_(inp), par_(par), cmd_(cmd), ang_(ang), verbose_(true)
{
    rad_.verbose(false);
    rad_.setupOneElectronIntegrals(par, cmd);
}

void Amplitudes::extract (std::string directory)
{
    par_.wait();

    if (verbose_)
    {
        std::cout << std::endl;
        if (cmd_.extract_extrapolate)
            std::cout << "Extrapolating T-matrices ";
        else
            std::cout << "Averaging T-matrices ";
        if (cmd_.extract_rho_begin > 0)
            std::cout << "from " << cmd_.extract_rho_begin << " ";
        if (cmd_.extract_rho > 0)
            std::cout << "to " << cmd_.extract_rho << " ";
        if (cmd_.extract_samples > 0)
            std::cout << "using " << cmd_.extract_samples << " samples";
        std::cout << std::endl;
    }

    for (unsigned Spin : inp_.Spin)
    {
        if (verbose_)
        {
            if (Spin == 0)
                std::cout << "\tSinglet" << std::endl;

            if (Spin == 1)
                std::cout << "\tTriplet" << std::endl;
        }

        for (unsigned ie = 0; ie < inp_.Etot.size(); ie++)
        {
            if (verbose_) std::cout << "\t\tEi = " << inp_.Etot[ie] << std::endl;

            // for all initial states
            for (auto instate  : inp_.instates)
            {
                // get initial quantum numbers
                int ni = std::get<0>(instate);
                int li = std::get<1>(instate);
                int mi = std::get<2>(instate);

                // is this solution allowed at all for the given angular basis?
                bool allowed = false;

                // -> find a valid combination of atomic and projectile angular momentum
                for (int l = std::abs(li - inp_.L); l <= li + inp_.L; l++)
                {
                    // does this combination conserve parity?
                    if ((inp_.L + li + l) % 2 != inp_.Pi)
                        continue;

                    // does this combination have valid 'mi' for this partial wave?
                    if (special::ClebschGordan(li,mi,l,0,inp_.L,mi) != 0)
                        allowed = true;
                }

                // skip forbidden wave functions
                if (not allowed)
                    continue;

                // check existence of the solution; take into account distributed calculations
                reader_ = SolutionIO (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Etot[ie], ang_, std::vector<std::pair<int,int>>(), directory + "/psi");
                BlockArray<Complex> solution (ang_.size(), true, "sol");
                std::size_t valid_blocks = 0;
                for (unsigned ill = 0; ill < ang_.size(); ill++) if (par_.isMyGroupWork(ill) and (not cmd_.shared_scratch or par_.IamGroupMaster()))
                {
                    if (reader_.load(solution, ill))
                        valid_blocks++;
                }
                par_.syncsum(&valid_blocks, 1);

                if (valid_blocks != ang_.size())
                {
                    if (verbose_)
                        std::cout << "\t\t\tSolution files for L = " << inp_.L << ", Pi = " << inp_.Pi << ", (ni,li,mi) = (" << ni << "," << li << "," << mi << ") not found." << std::endl;

                    continue;
                }

                // extract amplitudes to all final states
                for (auto outstate : inp_.outstates)
                {
                    int nf = std::get<0>(outstate);
                    int lf = std::get<1>(outstate);

                    // skip angular forbidden states
                    if (not allowed)
                        continue;

                    if (nf > 0)
                    {
                        //
                        // Discrete transition
                        //

                        if (verbose_)
                            std::cout << format("\t\t\texc: (%d,%d,%d) -> (%d,%d,*) ",ni, li, mi, nf, lf) << std::endl;

                        // compute radial integrals
                        for (int mf = -lf; mf <= lf; mf++)
                        {
                            // transition
                            Transition transition = { ni, li, mi, nf, lf, mf };

                            // extract radial part of the amplitude
                            computeLambda_(transition, solution, ie, Spin);
                        }
                    }
                    else if (inp_.Etot[ie] > 0)
                    {
                        //
                        // Ionization
                        //

                        if (verbose_)
                            std::cout << format("\t\t\tion: (%d,%d,%d) -> ion ",ni, li, mi) /* << std::endl */;

                        // transition
                        Transition transition = { ni, li, mi, 0, 0, 0 };

                        // compute Ξ
                        computeXi_(transition, solution, ie, Spin);
                    }
                }
            }
        }
    }

    // update T-matrices and cross sections for discrete transitions
    for (auto lambda : Lambda_Slp)
    {
        computeTmat_(lambda.first);
        computeSigma_(lambda.first);
    }

    // update cross sections for ionization
    for (auto xi : Xi_Sl1l2)
    {
        computeSigmaIon_(xi.first);
    }
}

void Amplitudes::writeSQL_files (std::string directory)
{
    // let the file be written by the master process, and by the group-masters in multi-group case
    if (not par_.IamMaster() and (cmd_.shared_scratch or not par_.IamGroupMaster()))
        return;

    // compose output filename
    std::ostringstream ossfile;
#ifdef WITH_BOINC
    ossfile << "tmat.sql";
#else
    ossfile << "tmat-L" << inp_.L << "-Pi" << inp_.Pi << ".sql";
#endif

    // Create SQL batch file
    std::ofstream fsql (directory + "/" + ossfile.str().c_str());

    // set exponential format for floating point output
    fsql.setf(std::ios_base::scientific);

    // write header
    fsql << logo("--");
    fsql << "-- File generated on " << current_time();
    fsql << "--" << std::endl;
    fsql << "-- Partial T-matrices for use in the database interface program \"hex-db\"." << std::endl;
    fsql << "-- Use for example:" << std::endl;
    fsql << "--    > hex-db --new --database hex.db --import " << ossfile.str().c_str() << " --update" << std::endl;
    fsql << "--" << std::endl;
    fsql << "BEGIN TRANSACTION;" << std::endl;

    // for all discrete transitions data
    for (auto Tmat : Tmat_Slp)
    {
        // get transition
        Transition const & T = Tmat.first;

        // for all angular momenta (partial waves)
        for (int ell = 0; ell <= inp_.maxell; ell++)
        {
            // get T-matrices
            cArray const & T_S0 = Tmat.second[ell].first;
            cArray const & T_S1 = Tmat.second[ell].second;

            // write energies
            for (unsigned i = 0; i < inp_.Etot.size(); i++)
            {
                // write singlet value (S = 0)
                if (Complex_finite(T_S0[i]) and T_S0[i] != 0.0_r)
                {
                    fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                        << T.ni << "," << T.li << "," << T.mi << ","
                        << T.nf << "," << T.lf << "," << T.mf << ","
                        << inp_.L  << "," << 0 << ","
                        << inp_.Etot[i] + 1. / (T.ni * T.ni) << "," << ell << "," 
                        << T_S0[i].real() << "," << T_S0[i].imag() << ");" << std::endl;
                }

                // write triplet value (S = 1)
                if (Complex_finite(T_S1[i]) and T_S1[i] != 0.0_r)
                {
                    fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                        << T.ni << "," << T.li << "," << T.mi << ","
                        << T.nf << "," << T.lf << "," << T.mf << ","
                        << inp_.L  << "," << 1 << ","
                        << inp_.Etot[i] + 1. / (T.ni * T.ni) << "," << ell << "," 
                        << T_S1[i].real() << "," << T_S1[i].imag() << ");" << std::endl;
                }
            }
        }
    }

    // for all ionizations
    for (auto xi : Xi_Sl1l2)
    {
        // get transition
        Transition const & T = xi.first;
        std::vector<std::pair<cArray,cArray>> const & data = xi.second;

        // for all energies and angular momenta
        for (unsigned ill = 0; ill < ang_.size(); ill++) //??? or triangular
        for (std::size_t ie = 0; ie < inp_.Etot.size(); ie++)
        {
            // get Chebyshev expansion coefficients
            cArray const & Xi_S0 = data[ill * inp_.Etot.size() + ie].first;
            cArray const & Xi_S1 = data[ill * inp_.Etot.size() + ie].second;

            // save singlet data as BLOBs
            fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                 << T.ni << "," << T.li << "," << T.mi << ","
                 << inp_.L << ", 0, " << inp_.Etot[ie] + 1. / (T.ni * T.ni) << ", "
                 << ang_[ill].first << ", " << ang_[ill].second << ", "
                 << Xi_S0.toBlob().c_str() << ");" << std::endl;

            // save triplet data as BLOBs
            fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                 << T.ni << "," << T.li << "," << T.mi << ","
                 << inp_.L << ", 1, " << inp_.Etot[ie] + 1. / (T.ni * T.ni) << ", "
                 << ang_[ill].first << ", " << ang_[ill].second << ", "
                 << Xi_S1.toBlob().c_str() << ");" << std::endl;
        }
    }

    // finish writing
    fsql << "COMMIT;" << std::endl;
    fsql.close();
}

void Amplitudes::writeICS_files (std::string directory)
{
    // let the file be written by the master process, and by the group-masters in multi-group case
    if (not par_.IamMaster() and (cmd_.shared_scratch or not par_.IamGroupMaster()))
        return;

    // open files
    std::ofstream fS0 (directory + "/" + format("ics-L%d-S0-Pi%d.dat", inp_.L, inp_.Pi));
    std::ofstream fS1 (directory + "/" + format("ics-L%d-S1-Pi%d.dat", inp_.L, inp_.Pi));

    // write singlet file header
    fS0 << logo("#");
    fS0 << "# File generated on " << current_time() << "#" << std::endl;
    fS0 << "# Singlet partial cross sections." << std::endl  << "#" << std::endl;

    // write triplet file header
    fS1 << logo("#");
    fS1 << "# File generated on " << current_time() << "#" << std::endl;
    fS1 << "# Triplet partial cross sections." << std::endl << "#" << std::endl;

    // print column headers
    fS0 << std::left << std::setw(15) << "# E[Ry]";
    fS1 << std::left << std::setw(15) << "# E[Ry]";
    for (auto data : sigma_S)
    {
        // get transition
        Transition const & T = data.first;

        // write transition
        std::string header = format
        (
            "%s-%s",
            Hydrogen::stateName(T.ni,T.li,T.mi).c_str(),
            Hydrogen::stateName(T.nf,T.lf,T.mf).c_str()
        );
        fS0 << std::setw(15) << header;
        fS1 << std::setw(15) << header;
    }
    fS0 << std::endl << std::setw(15) << "# ----------";
    fS1 << std::endl << std::setw(15) << "# ----------";
    for (auto data : sigma_S)
    {
        fS0 << std::setw(15) << "----------";
        fS1 << std::setw(15) << "----------";
    }
    fS0 << std::endl;
    fS1 << std::endl;

    // print data (cross sections)
    for (unsigned ie = 0; ie < inp_.Etot.size(); ie++)
    {
        fS0 << std::setw(15) << inp_.Etot[ie];
        fS1 << std::setw(15) << inp_.Etot[ie];

        for (auto data : sigma_S)
        {
            // get singlet and triplet data
            rArray const & sigma_S0 = data.second.first;
            rArray const & sigma_S1 = data.second.second;

            fS0 << std::setw(15) << (std::isfinite(sigma_S0[ie]) ? sigma_S0[ie] : 0.);
            fS1 << std::setw(15) << (std::isfinite(sigma_S1[ie]) ? sigma_S1[ie] : 0.);
        }
        fS0 << std::endl;
        fS1 << std::endl;
    }

    // finish writing
    fS0.close();
    fS1.close();
}

cArray Amplitudes::readAtomPseudoState (int l, int ichan) const
{
    std::string filename = format("Hl-1-%d-%.4x.hdf", l, rad_.bspline_x().hash());

    HDFFile datafile (filename, HDFFile::readonly);
    if (not datafile.valid())
        HexException("File %s with the one-electron eigenstates was not found. Run the solver again to regenerate it.", filename.c_str());

    int N = datafile.size("Dl") / 2;
    if (N != rad_.bspline_x().Nspline())
        HexException("File %s is not compatible with the current basis. You should delete it.", filename.c_str());
    if (ichan >= N)
        HexException("File %s contains only %d (<= %d) channels.", filename.c_str(), N, ichan);

    cArray energies (N);
    if (not datafile.read("Dl", &energies[0], N))
        HexException("File %s does not contain the requested dataset \"Dl\".", filename.c_str());

    iArray indices (N);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j){ return energies[i].real() < energies[j].real(); });

    unsigned Nspline_inner = rad_.bspline_x().Nspline();

    cArray data (Nspline_inner);
    if (not datafile.read("Cl", &data[0], Nspline_inner, indices[ichan] * Nspline_inner))
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

void Amplitudes::computeLambda_ (Amplitudes::Transition T, BlockArray<Complex> & solution, int ie, int Spin)
{
    // final projectile energy
    Real Ef = inp_.Etot[ie] + 1.0_r/(T.nf*T.nf) + (T.mf-T.mi) * inp_.B;

    // final projectile momentum
    Real kf = (Ef >= 0 ? std::sqrt(Ef) : special::constant::Nan);

    // shorthands
    unsigned Nenergy = inp_.Etot.size();               // energy count
    int order   = inp_.order;                   // B-spline order
    int Nspline_inner = bspline_inner_.Nspline(); // B-spline count (inner basis)
    int Nspline_full  = bspline_full_ .Nspline(); // B-spline count (combined basis)
    int Nspline_outer = Nspline_full - Nspline_inner; // B-spline count (outer basis)

    // check that memory for this transition is allocated
    if (Lambda_Slp.find(T) == Lambda_Slp.end())
    {
        // create new entry for this transition
        Lambda_Slp[T] = std::vector<std::pair<cArray,cArray>>(inp_.maxell + 1);

        // allocate arrays for all partial waves
        for (int ell = 0; ell <= inp_.maxell; ell++)
            Lambda_Slp[T][ell].first = Lambda_Slp[T][ell].second = cArray(Nenergy);
    }

    // skip impact energies with undefined outgoing momentum
    if (not std::isfinite(kf) or kf == 0.)
        return;

    // skip angular-forbidden transitions
    if (T.li > inp_.maxell or T.lf > inp_.maxell)
        return;

    // read the appropriate projectile channel function for the final state (nf,lf,*)
    cArray Xp, Sp;
    if (cmd_.analytic_eigenstates)
    {
        Sp = RadialIntegrals::overlapP(bspline_inner_, rad_.gaussleg(), inp_.Za, T.nf, T.lf);
    }
    else
    {
        Xp = readAtomPseudoState(T.lf, T.nf - T.lf - 1);
        Sp = rad_.S_x().dot(Xp);
    }

    // The extracted T-matrix oscillates and slowly radially converges.
    // If we are far enough and only oscillations are left, we can average several uniformly spaced
    // extractions and get rid of the oscillations. Otherwise we need to extrapolate also
    // the trend of the T-matrix.

    Real wavelength = special::constant::two_pi / kf;
    Real Rb     = (cmd_.extract_rho       > 0) ? cmd_.extract_rho       : bspline_full_.R2();
    Real Ra     = (cmd_.extract_rho_begin > 0) ? cmd_.extract_rho_begin : Rb - wavelength; Ra = std::max(bspline_full_.R2() / 2, Ra);
    int samples = (cmd_.extract_samples   > 0) ? cmd_.extract_samples   : 10;

    rArray grid;
    cArrays singlet_lambda, triplet_lambda;
    for (int ell = 0; ell <= inp_.maxell; ell++)
    {
        // resize arrays
        singlet_lambda.push_back(cArray(samples));
        triplet_lambda.push_back(cArray(samples));
    }

    if (Ra > Rb)
        HexException("Wrong order of radial extraction bounds, %g > %g.", Ra, Rb);

    if (Rb > bspline_full_.R2())
        HexException("Extraction radius too far, %g > %g.", Rb, bspline_full_.R2());

    // evaluate radial part for all evaluation radii
    for (int i = 0; i < samples; i++)
    {
        // this is the evaluation point
        Real eval_r = (samples == 1 ? Rb : Ra + (i + 1) * (Rb - Ra) / (samples + 1));
        grid.push_back(eval_r);

        // determine knot
        int eval_knot = bspline_full_.knot(eval_r);

        // evaluate j and dj at far radius for all angular momenta up to maxell
        // FIXME : Coulomb for charge (inp_->Za - 1)
        cArray j_R0 = special::ric_jv(inp_.maxell, kf * eval_r);
        cArray dj_R0 = special::dric_jv(inp_.maxell, kf * eval_r) * kf;

        // evaluate B-splines and their derivatives at evaluation radius
        cArray Bspline_R0 (Nspline_full), Dspline_R0 (Nspline_full);
        for (int ispline = 0; ispline < Nspline_full; ispline++)
        {
            // evaluate B-spline
            Bspline_R0[ispline] = bspline_full_.bspline(ispline, eval_knot, order, eval_r);

            // evaluate B-spline derivative
            Dspline_R0[ispline] = bspline_full_.dspline(ispline, eval_knot, order, eval_r);
        }

        // evaluate Wronskians
        cArrays Wj (inp_.maxell + 1);
        for (int l = 0; l <= inp_.maxell; l++)
            Wj[l] = dj_R0[l] * Bspline_R0 - j_R0[l] * Dspline_R0;

        // compute radial factor
        # pragma omp parallel for schedule (dynamic,1) if (cmd_.parallel_extraction)
        for (unsigned ill = 0; ill < ang_.size(); ill++) if (par_.isMyGroupWork(ill) and par_.IamGroupMaster())
        {
            // skip blocks that do not contribute to (l1 = ) lf
            if (ang_[ill].first != T.lf)
                continue;

            // get projectile angular momentum
            int ell = ang_[ill].second;

            // evaluate radial integrals
            Complex lambda = 0;
            if (inp_.inner_only)
            {
                // change view to row-major dense matrix
                RowMatrixView<Complex> PsiSc
                (
                    Nspline_inner,  // rows
                    Nspline_inner,  // columns
                    solution[ill]   // data
                );

                // calculate radial integral
                lambda = (Sp | PsiSc | Wj[ell]);
            }
            else
            {

                // change view to row-major dense matrix
                cArrayView PsiScFf
                (
                    solution[ill],  // data
                    Nspline_inner * Nspline_inner + (reader_.channels()[ill].first + T.nf - T.lf - 1) * Nspline_outer, // offset
                    Nspline_outer   // elements
                );

                // calculate radial integral
                lambda = (PsiScFf | Wj[ell].slice(Nspline_inner, Nspline_full));
            }

            // update the stored value
            # pragma omp critical
            if (Spin == 0)
                singlet_lambda[ell][i] += lambda;
            else
                triplet_lambda[ell][i] += lambda;
        }
    }

    for (int ell = 0; ell <= inp_.maxell; ell++)
    {
        if (cmd_.extract_extrapolate)
        {
            // radial extrapolation
            Lambda_Slp[T][ell].first[ie] += inv_power_extrapolate(grid, singlet_lambda[ell]);
            Lambda_Slp[T][ell].second[ie] += inv_power_extrapolate(grid, triplet_lambda[ell]);
        }
        else
        {
            // plain averaging
            Lambda_Slp[T][ell].first[ie] += sum(singlet_lambda[ell]) / Real(samples);
            Lambda_Slp[T][ell].second[ie] += sum(triplet_lambda[ell]) / Real(samples);
        }
    }
}

void Amplitudes::computeTmat_ (Amplitudes::Transition T)
{
    // final projectile momenta
    // final projectile energies
    rArray Ef = inp_.Etot + 1.0_r/(T.nf*T.nf) + (T.mf-T.mi) * inp_.B;

    // final projectile momenta
    rArray inv_kf; for (Real ef : Ef) inv_kf.push_back(ef > 0 ? 1./std::sqrt(ef) : 0.);

    // allocate memory
    if (Tmat_Slp.find(T) == Tmat_Slp.end())
    {
        Tmat_Slp[T] = std::vector<std::pair<cArray,cArray>>(inp_.maxell + 1);
        for (auto vec : Tmat_Slp[T])
            vec.first = vec.second = cArray(inp_.Etot.size());
    }

    // for all radial integrals (indexed by angular momenta)
    for (int ell = 0; ell <= inp_.maxell; ell++)
    {
        // get radial integrals
        cArray & rad_S0 = Lambda_Slp[T][ell].first;
        cArray & rad_S1 = Lambda_Slp[T][ell].second;

        // synchronize data
        par_.mastersum(rad_S0.data(), rad_S0.size(), 0);
        par_.mastersum(rad_S1.data(), rad_S1.size(), 0);
        par_.bcast(0, rad_S0);
        par_.bcast(0, rad_S1);

        // symmetry factor
        Real sf = (inp_.Zp > 0 ? 1.0_r : special::constant::sqrt_half);

        // compute T-matrices
        Tmat_Slp[T][ell].first = rad_S0 * 4.0_r * special::constant::pi * inv_kf * std::pow(Complex(0.,1.), -ell)
                    * (Real)special::ClebschGordan(T.lf, T.mf, ell, T.mi - T.mf, inp_.L, T.mi) * sf;
        Tmat_Slp[T][ell].second = rad_S1 * 4.0_r * special::constant::pi * inv_kf * std::pow(Complex(0.,1.), -ell)
                    * (Real)special::ClebschGordan(T.lf, T.mf, ell, T.mi - T.mf, inp_.L, T.mi) * sf;
    }
}

void Amplitudes::computeSigma_ (Amplitudes::Transition T)
{
    // initial and final projectile energies
    rArray Ei = inp_.Etot + 1.0_r/(T.ni*T.ni);
    rArray Ef = inp_.Etot + 1.0_r/(T.nf*T.nf) + (T.mf-T.mi) * inp_.B;

    // final projectile momenta
    rArray inv_ki; for (Real ei : Ei) inv_ki.push_back(ei > 0 ? 1./std::sqrt(ei) : 0.);
    rArray     kf; for (Real ef : Ef)     kf.push_back(ef > 0 ?    std::sqrt(ef) : 0.);

    // allocate memory
    if (sigma_S.find(T) == sigma_S.end())
        sigma_S[T] = std::make_pair(rArray(inp_.Etot.size()),rArray(inp_.Etot.size()));

    // for all T-matrices (indexed by angular momenta)
    for (auto tmat : Tmat_Slp[T])
    {
        // get radial integrals
        cArray const & Tmat_S0 = tmat.first;
        cArray const & Tmat_S1 = tmat.second;

        // compute singlet contribution
        rArray Re_f0_ell = -realpart(Tmat_S0) / special::constant::two_pi;
        rArray Im_f0_ell = -imagpart(Tmat_S0) / special::constant::two_pi;
        sigma_S[T].first += 0.25 * kf * inv_ki * (Re_f0_ell * Re_f0_ell + Im_f0_ell * Im_f0_ell);

        // compute triplet contribution
        rArray Re_f1_ell = -realpart(Tmat_S1) / special::constant::two_pi;
        rArray Im_f1_ell = -imagpart(Tmat_S1) / special::constant::two_pi;
        sigma_S[T].second += 0.75 * kf * inv_ki * (Re_f1_ell * Re_f1_ell + Im_f1_ell * Im_f1_ell);
    }
}

Chebyshev<double,Complex> Amplitudes::fcheb (cArrayView const & PsiSc, Real kmax, int l1, int l2)
{
    // shorthands
    Complex const * const t = &(bspline_inner_.t(0));
    int Nspline = bspline_inner_.Nspline();
    int Nreknot = bspline_inner_.Nreknot();
    int order   = bspline_inner_.order();

    // determine evaluation radius
    Real rho = (cmd_.extract_rho > 0) ? cmd_.extract_rho : t[Nreknot-2].real();

    // debug output
    std::ofstream dbg ("debug.log");

    // we want to approximate the following function f_{ℓ₁ℓ₂}^{LS}(k₁,k₂)
    auto fLSl1l2k1k2 = [&](Real k1) -> Complex
    {
        if (k1 == 0 or k1*k1 >= kmax*kmax)
            return 0.;

        // compute momentum of the other electron
        Real k2 = std::sqrt(kmax*kmax - k1*k1);

        // Xi integrand
        auto integrand = [&](Real alpha) -> Complex
        {

            // precompute projectors
            Real cos_alpha = (alpha == special::constant::pi_half) ? 0. : std::cos(alpha);
            Real sin_alpha = std::sin(alpha);

            // precompute coordinates
            Real r1 = rho * cos_alpha;
            Real r2 = rho * sin_alpha;

            // evaluate Coulomb wave functions and derivatives
            double F1, F2, F1p, F2p;
            int err1 = special::coul_F(l1,k1,r1, F1,F1p);
            int err2 = special::coul_F(l2,k2,r2, F2,F2p);
            if (err1 != GSL_SUCCESS or err2 != GSL_SUCCESS)
            {
                std::cerr << "Errors while evaluating Coulomb function:" << std::endl;
                std::cerr << "\terr1 = " << err1 << std::endl;
                std::cerr << "\terr2 = " << err2 << std::endl;
                std::exit(EXIT_FAILURE);
            }

            Real F1F2 = F1 * F2;
            Real ddrho_F1F2 = 0.;

            if (cos_alpha != 0.)
                ddrho_F1F2 += k1*F1p*cos_alpha*F2;
            if (sin_alpha != 0.)
                ddrho_F1F2 += k2*F1*F2p*sin_alpha;

            // get B-spline knots
            int iknot1 = bspline_inner_.knot(r1);
            int iknot2 = bspline_inner_.knot(r2);

            // auxiliary variables
            cArray B1 (Nspline), dB1 (Nspline), B2 (Nspline), dB2 (Nspline);

            // evaluate the B-splines
            for (int ispline1 = std::max(0,iknot1-order); ispline1 <= iknot1; ispline1++)
            {
                B1[ispline1]  = bspline_inner_.bspline(ispline1,iknot1,order,r1);
                dB1[ispline1] = bspline_inner_.dspline(ispline1,iknot1,order,r1);
            }
            for (int ispline2 = std::max(0,iknot2-order); ispline2 <= iknot2; ispline2++)
            {
                B2[ispline2]  = bspline_inner_.bspline(ispline2,iknot2,order,r2);
                dB2[ispline2] = bspline_inner_.dspline(ispline2,iknot2,order,r2);
            }

            // evaluate the solution
            Complex Psi = 0., ddr1_Psi = 0., ddr2_Psi = 0., ddrho_Psi = 0.;
            for (int ispline1 = std::max(0,iknot1-order); ispline1 <= iknot1; ispline1++)
            for (int ispline2 = std::max(0,iknot2-order); ispline2 <= iknot2; ispline2++)
            {
                int idx = ispline1 * Nspline + ispline2;

                Psi      += PsiSc[idx] *  B1[ispline1] *  B2[ispline2];
                ddr1_Psi += PsiSc[idx] * dB1[ispline1] *  B2[ispline2];
                ddr2_Psi += PsiSc[idx] *  B1[ispline1] * dB2[ispline2];
            }

            if (cos_alpha != 0.)
                ddrho_Psi += ddr1_Psi * cos_alpha;
            if (sin_alpha != 0.)
                ddrho_Psi += ddr2_Psi * sin_alpha;

            /// DEBUG
            if (not std::isfinite(F1F2))
                std::cerr << "F1F2 = " << F1F2 << std::endl;
            if (not std::isfinite(std::abs(ddrho_Psi)))
                std::cerr << "ddrho_Psi = " << ddrho_Psi << std::endl;
            if (not std::isfinite(std::abs(Psi)))
                std::cerr << "Psi = " << Psi << std::endl;
            if (not std::isfinite(ddrho_F1F2))
                std::cerr << "ddrho_F1F2 = " << ddrho_F1F2 << std::endl;
            ///

            // evaluate the integrand
            return F1F2*ddrho_Psi - Psi*ddrho_F1F2;
        };

        // integrator
        ClenshawCurtis<decltype(integrand),Complex> Q(integrand);
        Q.setEps(1e-6);
        Complex res = 2.0_r * rho * Q.integrate(0., special::constant::pi_half) / special::constant::sqrt_pi;

        dbg << "CB " << k1 << " " << res.real() << " " << res.imag() << std::endl;

        return res;
    };

    // Chebyshev approximation
    Chebyshev<double,Complex> CB;

    // avoid calculation when the extraction radius is too far
    if (rho > bspline_inner_.R2())
    {
        std::cout << "Warning: Extraction radius rho = " << rho << " is too far; the atomic real grid ends at R0 = " << bspline_inner_.R2() << std::endl;
    }
    else   
    {
        // convergence loop
        for (int N = 4; ; N *= 2)
        {
            // build the approximation
            CB.generate(fLSl1l2k1k2, N, 0., kmax);

            // check tail
            if (CB.tail(1e-5) != N)
                break;

            // limit subdivision
            if (N > 32768)
                HexException("ERROR: Non-convergent Chebyshev expansion.");
        }
    }

    return CB;
}

void Amplitudes::computeXi_ (Amplitudes::Transition T, BlockArray<Complex> & solution, int ie, int Spin)
{
    // maximal available momentum
    double kmax = std::sqrt(inp_.Etot[ie]);

    // check that memory for this transition is allocated; allocate if not
    if (Xi_Sl1l2.find(T) == Xi_Sl1l2.end())
        Xi_Sl1l2[T] = std::vector<std::pair<cArray,cArray>>(ang_.size() * inp_.Etot.size());

    if (verbose_) std::cout << ": ";

    // for all angular states
    for (unsigned ill = 0; ill < ang_.size(); ill++)
    {
        int l1 = ang_[ill].first;
        int l2 = ang_[ill].second;

        if (verbose_) std::cout << "(" << l1 << "," << l2 << ") " << std::flush;

        // load solution block
        if (not solution.inmemory())
            solution.hdfload(ill);

        // compute new ionization amplitude
        Chebyshev<double,Complex> CB = fcheb(solution[ill], kmax, l1, l2);
        if (Spin == 0)
            Xi_Sl1l2[T][ill * inp_.Etot.size() + ie].first = CB.coeffs();
        else
            Xi_Sl1l2[T][ill * inp_.Etot.size() + ie].second = CB.coeffs();

        // unload solution block
        if (not solution.inmemory())
            solution[ill].drop();
    }

    if (verbose_) std::cout << "ok" << std::endl;
}

void Amplitudes::computeSigmaIon_ (Amplitudes::Transition T)
{
    // number of energies
    unsigned Nenergy = inp_.Etot.size();

    // initial momentum
    rArray ki = sqrt(inp_.Etot + 1.0_r/(T.ni*T.ni));

    // allocate memory for cross sections
    if (sigma_S.find(T) == sigma_S.end())
        sigma_S[T] = std::make_pair(rArray(Nenergy),rArray(Nenergy));

    // for all energies and angular blocks
    for (std::size_t ie = 0; ie < Nenergy; ie++)
    for (unsigned ill = 0; ill < ang_.size(); ill++)
    {
        // maximal available momentum
        double kmax = std::sqrt(inp_.Etot[ie]);

        // Chebyshev expansion coefficients
        Chebyshev<double,Complex> CB;

        // integrand |f|²
        int tail; int n;
        auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(kmax * std::sin(beta), tail)); };

        // integrator
        ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);

        // integrate singlet
        if (not Xi_Sl1l2[T][ill * Nenergy + ie].first.empty())
        {
            CB = Chebyshev<double,Complex>(Xi_Sl1l2[T][ill * Nenergy + ie].first, 0., kmax); tail = CB.tail(1e-10); 
            sigma_S[T].first[ie] += integrator.integrate(0, special::constant::pi_quart, &n) / ki[ie];
        }

        // integrate triplet
        if (not Xi_Sl1l2[T][ill * Nenergy + ie].second.empty())
        {
            CB = Chebyshev<double,Complex>(Xi_Sl1l2[T][ill * Nenergy + ie].second, 0., kmax); tail = CB.tail(1e-10); 
            sigma_S[T].second[ie] += integrator.integrate(0, special::constant::pi_quart, &n) / ki[ie];
        }
    }
}
