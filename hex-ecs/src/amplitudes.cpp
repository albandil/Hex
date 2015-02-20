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

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "chebyshev.h"
#include "clenshawcurtis.h"
#include "hydrogen.h"
#include "matrix.h"
#include "radial.h"
#include "special.h"
#include "version.h"

std::string current_time () 
{
    std::time_t result = std::time(nullptr);
    return std::asctime(std::localtime(&result));
}

void Amplitudes::extract ()
{
    // radial integrals
    RadialIntegrals rad (bspline_);
    
    std::cout << std::endl << "Extracting T-matrices" << std::endl;
    
    for (unsigned Spin = 0; Spin <= 1; Spin++)
    {
        if (Spin == 0)
            std::cout << "\tSinglet" << std::endl;
        
        if (Spin == 1)
            std::cout << "\tTriplet" << std::endl;
        
        for (unsigned ie = 0; ie < inp_.Ei.size(); ie++)
        {
            std::cout << "\t\tEi = " << inp_.Ei[ie] << std::endl;
            
            // for all initial states
            for (auto instate  : inp_.instates)
            {
                // get initial quantum numbers
                int ni = std::get<0>(instate);
                int li = std::get<1>(instate);
                int mi = std::get<2>(instate);
                
                // load the solution
                SolutionIO reader (inp_.L, Spin, inp_.Pi, ni, li, mi, inp_.Ei[ie], ang_, bspline_.Nspline());
                cArray solution;
                if (not reader.load(solution))
                {
                    std::cout << "\t\t\tFile \"" << reader.name() << "\" not found." << std::endl;
                    continue;
                }
                
                // extract amplitudes to all final states
                for (auto outstate : inp_.outstates)
                {
                    int nf = std::get<0>(outstate);
                    int lf = std::get<1>(outstate);
                    
                    // check if the right hand side will be zero for this instate
                    bool allowed = false;
                    for (int l = std::abs(li - inp_.L); l <= li + inp_.L; l++)
                    {
                        // does this combination conserve parity?
                        if ((inp_.L + li + l) % 2 != inp_.Pi)
                            continue;
                        
                        // does this combination have valid 'mi' for this partial wave?
                        if (special::ClebschGordan(li,mi,l,0,inp_.L,mi) != 0)
                            allowed = true;
                    }
                    
                    // skip angular forbidden states
                    if (not allowed)
                        continue;
                    
                    if (nf > 0)
                    {
                        //
                        // Discrete transition
                        //
                        
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
                    else
                    {
                        //
                        // Ionization
                        //
                        
                        std::cout << format("\tion: (%d,%d,%d) -> ion ",ni, li, mi) << std::endl;
                        
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

void Amplitudes::writeSQL_files ()
{
    // compose output filename
    std::ostringstream ossfile;
    if (par_.active())
    {
        ossfile << "tmat-n" << inp_.ni << "-L" << inp_.L << "-Pi" << inp_.Pi << "-(" << par_.iproc() << ").sql";
    }
    else
    {
        ossfile << "tmat-n" << inp_.ni << "-L" << inp_.L << "-Pi" << inp_.Pi << ".sql";
    }
    
    // Create SQL batch file
    std::ofstream fsql(ossfile.str().c_str());
        
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
            for (unsigned i = 0; i < inp_.Ei.size(); i++)
            {
                // write singlet value (S = 0)
                if (Complex_finite(T_S0[i]) and T_S0[i] != 0.)
                {
                    fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                        << T.ni << "," << T.li << "," << T.mi << ","
                        << T.nf << "," << T.lf << "," << T.mf << ","
                        << inp_.L  << "," << 0 << ","
                        << inp_.Ei[i] << "," << ell << "," 
                        << T_S0[i].real() << "," << T_S0[i].imag() << ","
                        << "0,0);" << std::endl;
                }
                
                // write triplet value (S = 1)
                if (Complex_finite(T_S1[i]) and T_S1[i] != 0.)
                {
                    fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                        << T.ni << "," << T.li << "," << T.mi << ","
                        << T.nf << "," << T.lf << "," << T.mf << ","
                        << inp_.L  << "," << 1 << ","
                        << inp_.Ei[i] << "," << ell << "," 
                        << T_S1[i].real() << "," << T_S1[i].imag() << ","
                        << "0,0);" << std::endl;
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
        for (std::size_t ie = 0; ie < inp_.Ei.size(); ie++)
        {
            // get Chebyshev expansion coefficients
            cArray const & Xi_S0 = data[ill * inp_.Ei.size() + ie].first;
            cArray const & Xi_S1 = data[ill * inp_.Ei.size() + ie].second;
            
            // save singlet data as BLOBs
            fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                 << inp_.ni << "," << T.li << "," << T.mi << ","
                 << inp_.L << ", 0, " << inp_.Ei[ie] << ", "
                 << ang_[ill].first << ", " << ang_[ill].second << ", "
                 << Xi_S0.toBlob().c_str() << ");" << std::endl;
            
            // save triplet data as BLOBs
            fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
                 << inp_.ni << "," << T.li << "," << T.mi << ","
                 << inp_.L << ", 1, " << inp_.Ei[ie] << ", "
                 << ang_[ill].first << ", " << ang_[ill].second << ", "
                 << Xi_S1.toBlob().c_str() << ");" << std::endl;
        }
    }
    
    // finish writing
    fsql << "COMMIT;" << std::endl;
    fsql.close();
}

void Amplitudes::writeICS_files ()
{
    // open files
    std::ofstream fS0 (format("ics-n%d-L%d-S0-Pi%d.dat", inp_.ni, inp_.L, inp_.Pi));
    std::ofstream fS1 (format("ics-n%d-L%d-S1-Pi%d.dat", inp_.ni, inp_.L, inp_.Pi));
    
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
    for (unsigned ie = 0; ie < inp_.Ei.size(); ie++)
    {
        fS0 << std::setw(15) << inp_.Ei[ie];
        fS1 << std::setw(15) << inp_.Ei[ie];
        
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

void Amplitudes::computeLambda_ (Amplitudes::Transition T, const cArrayView solution, int ie, int Spin)
{
    // final projectile momenta
    rArray kf = sqrt(inp_.Ei - 1./(inp_.ni*inp_.ni) + 1./(T.nf*T.nf) + (T.mf-T.mi) * inp_.B);
    
    // shorthands
    unsigned Nenergy = kf.size();                // energy count
    Complex const * const t = &(bspline_.t(0));   // B-spline knots
    int order   = bspline_.order();               // B-spline order
    int Nspline = bspline_.Nspline();             // B-spline count
    int Nreknot = bspline_.Nreknot();             // number of real knots
    
    // compute final hydrogen orbital overlaps with B-spline basis
    cArray Pf_overlaps = rad_.overlapP(T.nf, T.lf, weightEndDamp(bspline_));
    
    // check that memory for this transition is allocated
    if (Lambda_Slp.find(T) == Lambda_Slp.end())
    {
        // create new entry for this transition
        Lambda_Slp[T] = std::vector<std::pair<cArray,cArray>>(inp_.maxell + 1);
        
        // allocate arrays for all partial waves
        for (int ell = 0; ell <= inp_.maxell; ell++)
            Lambda_Slp[T][ell].first = Lambda_Slp[T][ell].second = cArray(Nenergy);
    }
    
    // The cross section oscillates, so we will do some averaging
    // As recommended by Bartlett, we will compute several amplitudes
    // separated by π/(n*kf[ie]) near the R₀ turning point.
    double wavelength = special::constant::pi / kf[ie];
    char const * HEX_RHO = std::getenv("HEX_RHO");
    char const * HEX_SAMPLES = std::getenv("HEX_SAMPLES");
    int samples = (HEX_SAMPLES == nullptr) ? 10 : std::atoi(HEX_SAMPLES);
    double R0 = (HEX_RHO == nullptr) ? t[Nreknot - 1].real() : std::atof(HEX_RHO);
    
    // skip impact energies with undefined outgoing momentum
    if (not std::isfinite(kf[ie]) or kf[ie] == 0.)
        return;
    
    // evaluate radial part for all evaluation radii
    # pragma omp parallel for
    for (int n = 1; n <= samples; n++)
    {
        // this is the evaluation point
        double eval_r = R0 - wavelength * n / samples;
        
        // determine knot
        int eval_knot = bspline_.knot(eval_r);
        
        // evaluate j and dj at far radius for all angular momenta up to maxell
        cArray j_R0 = special::ric_jv(inp_.maxell, kf[ie] * eval_r);
        cArray dj_R0 = special::dric_jv(inp_.maxell, kf[ie] * eval_r) * kf[ie];
        
        // evaluate B-splines and their derivatives at evaluation radius
        CooMatrix Bspline_R0(Nspline, 1), Dspline_R0(Nspline, 1);
        for (int ispline = 0; ispline < Nspline; ispline++)
        {
            Complex val;
            
            // evaluate B-spline
            if ((val = bspline_.bspline(ispline, eval_knot, order, eval_r)) != 0.)
                Bspline_R0.add(ispline, 0, val);
            
            // evaluate B-spline derivative
            if ((val = bspline_.dspline(ispline, eval_knot, order, eval_r)) != 0.)
                Dspline_R0.add(ispline, 0, val);
        }
        
        // evaluate Wronskians
        CooMatrix Wj[inp_.maxell + 1];
        for (int l = 0; l <= inp_.maxell; l++)
            Wj[l] = dj_R0[l] * Bspline_R0 - j_R0[l] * Dspline_R0;
        
        // we need "P_overlaps" to have a 'dot' method
        CooMatrix Sp (Nspline, 1, Pf_overlaps.begin());
        
        // compute radial factor
        for (unsigned ill = 0; ill < ang_.size(); ill++)
        {
            // skip blocks that do not contribute to (l1 = ) lf
            if (ang_[ill].first != T.lf)
                continue;
            
            // get angular momentum
            int ell = ang_[ill].second;
            
            // get correct solution (for this ang. mom.)
            cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
            
            // update value in the storage
            Complex lambda = Sp.transpose().dot(PsiSc).dot(Wj[ell].todense()).todense()[0] / double(samples);
            # pragma omp critical
            if (Spin == 0)
                (Lambda_Slp[T][ell].first)[ie] += lambda;
            else
                (Lambda_Slp[T][ell].second)[ie] += lambda;
        }
    }
}

void Amplitudes::computeTmat_ (Amplitudes::Transition T)
{
    // final projectile momenta
    rArray kf = sqrt(inp_.Ei - 1./(inp_.ni*inp_.ni) + 1./(T.nf*T.nf) + (T.mf-T.mi) * inp_.B);
    
    // allocate memory
    if (Tmat_Slp.find(T) == Tmat_Slp.end())
    {
        Tmat_Slp[T] = std::vector<std::pair<cArray,cArray>>(inp_.maxell + 1);
        for (auto vec : Tmat_Slp[T])
            vec.first = vec.second = cArray(inp_.Ei.size());
    }
    
    // for all radial integrals (indexed by angular momenta)
    for (int ell = 0; ell <= inp_.maxell; ell++)
    {
        // get radial integrals
        cArray const & rad_S0 = Lambda_Slp[T][ell].first;
        cArray const & rad_S1 = Lambda_Slp[T][ell].second;
        
        // compute T-matrices
        Tmat_Slp[T][ell].first = rad_S0 * 4. * special::constant::pi / kf * std::pow(Complex(0.,1.), -ell)
                    * special::ClebschGordan(T.lf, T.mf, ell, T.mi - T.mf, inp_.L, T.mi) * special::constant::sqrt_half;
        Tmat_Slp[T][ell].second = rad_S1 * 4. * special::constant::pi / kf * std::pow(Complex(0.,1.), -ell)
                    * special::ClebschGordan(T.lf, T.mf, ell, T.mi - T.mf, inp_.L, T.mi) * special::constant::sqrt_half;
    }
}

void Amplitudes::computeSigma_ (Amplitudes::Transition T)
{
    // final projectile momenta
    rArray kf = sqrt(inp_.Ei - 1./(inp_.ni*inp_.ni) + 1./(T.nf*T.nf) + (T.mf-T.mi) * inp_.B);
    
    // allocate memory
    if (sigma_S.find(T) == sigma_S.end())
        sigma_S[T] = std::make_pair(rArray(inp_.Ei.size()),rArray(inp_.Ei.size()));
    
    // for all T-matrices (indexed by angular momenta)
    for (auto tmat : Tmat_Slp[T])
    {
        // get radial integrals
        cArray const & Tmat_S0 = tmat.first;
        cArray const & Tmat_S1 = tmat.second;
        
        // compute singlet contribution
        rArray Re_f0_ell = -realpart(Tmat_S0) / special::constant::two_pi;
        rArray Im_f0_ell = -imagpart(Tmat_S0) / special::constant::two_pi;
        sigma_S[T].first += 0.25 * kf / inp_.ki * (Re_f0_ell * Re_f0_ell + Im_f0_ell * Im_f0_ell);
        
        // compute triplet contribution
        rArray Re_f1_ell = -realpart(Tmat_S1) / special::constant::two_pi;
        rArray Im_f1_ell = -imagpart(Tmat_S1) / special::constant::two_pi;
        sigma_S[T].second += 0.75 * kf / inp_.ki * (Re_f1_ell * Re_f1_ell + Im_f1_ell * Im_f1_ell);
    }
}

Chebyshev<double,Complex> fcheb (Bspline const & bspline, cArrayView const & PsiSc, double kmax, int l1, int l2)
{
    // shorthands
    Complex const * const t = &(bspline.t(0));   // B-spline knots
    int Nspline = bspline.Nspline();             // number of real knots
    int Nreknot = bspline.Nreknot();             // number of real knots
    int order   = bspline.order();               // B-spline order
    
    // determine evaluation radius
    char const * HEX_RHO = std::getenv("HEX_RHO");
    double rho = (HEX_RHO == nullptr) ? t[Nreknot-2].real() : std::atof(HEX_RHO);
    
    // we want to approximate the following function f_{ℓ₁ℓ₂}^{LS}(k₁,k₂)
    auto fLSl1l2k1k2 = [&](double k1) -> Complex
    {
        if (k1 == 0 or k1*k1 >= kmax*kmax)
            return 0.;
        
        // compute momentum of the other electron
        double k2 = std::sqrt(kmax*kmax - k1*k1);
        
        // Xi integrand
        auto integrand = [&](double alpha) -> Complex
        {
            
            // precompute projectors
            double cos_alpha = (alpha == special::constant::pi_half) ? 0. : std::cos(alpha);
            double sin_alpha = std::sin(alpha);
            
            // precompute coordinates
            double r1 = rho * cos_alpha;
            double r2 = rho * sin_alpha;
            
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
            
            double F1F2 = F1 * F2;
            double ddrho_F1F2 = 0.;
            
            if (cos_alpha != 0.)
                ddrho_F1F2 += k1*F1p*cos_alpha*F2;
            if (sin_alpha != 0.)
                ddrho_F1F2 += k2*F1*F2p*sin_alpha;
            
            // get B-spline knots
            int iknot1 = bspline.knot(r1);
            int iknot2 = bspline.knot(r2);
            
            // auxiliary variables
            cArray B1(Nspline), dB1(Nspline), B2(Nspline), dB2(Nspline);
            
            // evaluate the B-splines
            for (int ispline1 = std::max(0,iknot1-order); ispline1 <= iknot1; ispline1++)
            {
                B1[ispline1]  = bspline.bspline(ispline1,iknot1,order,r1);
                dB1[ispline1] = bspline.dspline(ispline1,iknot1,order,r1);
            }
            for (int ispline2 = std::max(0,iknot2-order); ispline2 <= iknot2; ispline2++)
            {
                B2[ispline2]  = bspline.bspline(ispline2,iknot2,order,r2);
                dB2[ispline2] = bspline.dspline(ispline2,iknot2,order,r2);
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
        Complex res = 2. * rho * Q.integrate(0., special::constant::pi_half) / special::constant::sqrt_pi;
        
        return res;
    };
    
    // Chebyshev approximation
    Chebyshev<double,Complex> CB;
    
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
    
    return CB;
}

void Amplitudes::computeXi_ (Amplitudes::Transition T, const cArrayView solution, int ie, int Spin)
{
    // B-spline count
    int Nspline = bspline_.Nspline();
    
    // maximal available momentum
    double kmax = std::sqrt(inp_.Ei[ie] - 1./(T.ni*T.ni));
    
    // for all angular states ???: (triangle ℓ₂ ≤ ℓ₁)
    for (unsigned ill = 0; ill < ang_.size(); ill++)
    {
        int l1 = ang_[ill].first;
        int l2 = ang_[ill].second;
        
        std::cout << "(" << l1 << "," << l2 << ") " << std::flush;
        
        // create subset of the solution
        cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
        
        // compute new ionization amplitude
        Chebyshev<double,Complex> CB = fcheb(bspline_, PsiSc, kmax, l1, l2);
        if (Spin == 0)
            Xi_Sl1l2[T][ill * inp_.Ei.size() + ie].first = CB.coeffs();
        else
            Xi_Sl1l2[T][ill * inp_.Ei.size() + ie].second = CB.coeffs();
    }
    
    std::cout << std::endl;
}

void Amplitudes::computeSigmaIon_ (Amplitudes::Transition T)
{
    // number of energies
    unsigned Nenergy = inp_.Ei.size();
    
    // allocate memory for cross sections
    if (sigma_S.find(T) == sigma_S.end())
        sigma_S[T] = std::make_pair(rArray(Nenergy),rArray(Nenergy));
    
    // for all energies and angular blocks
    for (std::size_t ie = 0; ie < inp_.Ei.size(); ie++)
    for (unsigned ill = 0; ill < ang_.size(); ill++)
    {
        // maximal available momentum
        double kmax = std::sqrt(inp_.Ei[ie] - 1./(T.ni*T.ni));
        
        // Chebyshev expansion coefficients
        Chebyshev<double,Complex> CB;
        
        // integrand |f|²
        int tail; int n;
        auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(kmax * std::sin(beta), tail)); };
        
        // integrator
        ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
        
        // integrate singlet
        CB = Chebyshev<double,Complex>(Xi_Sl1l2[T][ill * inp_.Ei.size() + ie].first, 0., kmax); tail = CB.tail(1e-10); 
        sigma_S[T].first[ie] += integrator.integrate(0, special::constant::pi_quart, &n) / std::sqrt(inp_.Ei[ie]);
        
        // integrate triplet
        CB = Chebyshev<double,Complex>(Xi_Sl1l2[T][ill * inp_.Ei.size() + ie].second, 0., kmax); tail = CB.tail(1e-10); 
        sigma_S[T].second[ie] = integrator.integrate(0, special::constant::pi_quart, &n) / std::sqrt(inp_.Ei[ie]);
    }
}
