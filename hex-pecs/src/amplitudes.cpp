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

#include <gsl/gsl_interp.h>

#include "hex-hydrogen.h"
#include "hex-version.h"

#include "amplitudes.h"

Amplitudes::Amplitudes
(
    CommandLine const & cmd,
    InputFile const & inp,
    AngularBasis const & ang,
    RadialBasis const & rad
) : cmd_(cmd), inp_(inp), ang_(ang), rad_(rad)
{
    // resize T-matrix arrays and initialize to zero
    T_matrices.resize(cmd_.Epert.size());
    for (std::size_t e = 0; e < cmd_.Epert.size(); e++)
    {
        T_matrices[e].resize(inp_.istates.size());
        for (std::size_t i = 0; i < inp_.istates.size(); i++)
        {
            T_matrices[e][i].resize(inp_.fstates.size());
            for (std::size_t f = 0; f < inp_.istates.size(); f++)
            {
                T_matrices[e][i][f].resize(ang_.maxell() + 1);
            }
        }
    }
}

void Amplitudes::extract
(
    unsigned istate,
    unsigned iEpert,
    Complex const * psi
)
{
    gsl_interp * spline = gsl_interp_alloc(gsl_interp_cspline, rad_.rgrid.size());
    gsl_interp_accel * acc = gsl_interp_accel_alloc();
    
    // initial state quantum numbers
    //int ni = inp_.istates[istate].n;
    //int li = inp.fstates[istate].l; // - not used
    int mi = inp_.istates[istate].m;
    
    // energy perturbation
    double Epert = cmd_.Epert[iEpert];
    
    for (std::size_t fstate = 0; fstate < inp_.fstates.size(); fstate++)
    {
        // final state euantum numbers
        int nf = inp_.fstates[fstate].n;
        int lf = inp_.fstates[fstate].l;
        int mf = inp_.fstates[fstate].m;
        
        // final projectile momentum
        double kf = std::sqrt(inp_.Etot + Epert + 1./(nf*nf));
        
        // skip in-accessible states
        if (not std::isfinite(kf))
            continue;
        
        // evaluate the final hydrogen function in grid points
        rArray P (rad_.rgrid.size());
        for (std::size_t i = 0; i < rad_.rgrid.size(); i++)
            P[i] = Hydrogen::P(nf, lf, rad_.rgrid[i], inp_.Z);
        
        // for all contributing angular blocks
        for (std::size_t iblock = 0; iblock < ang_.size(); iblock++)
        if ((unsigned)inp_.fstates[fstate].l == ang_.state(iblock).first)
        {
            rArray jre (rad_.rgrid.size()), jim (rad_.rgrid.size()), Jre (rad_.rgrid.size()), Jim (rad_.rgrid.size());
            
            // extraction radius
            double r_eval = rad_.rgrid.back(1);
            
            // partial wave (projectile final angular momentum)
            unsigned ell = ang_.state(iblock).second;
            
            // calculate overlaps of all products and the final radial state
            for (std::size_t i = 0; i < rad_.rgrid.size(); i++)
            {
                // project column to radial function
                for (std::size_t j = 0; j < rad_.rgrid.size(); j++)
                {
                    Complex prod = P[j] * psi[(iblock * rad_.Npts + i) * rad_.Npts + j];
                    jre[j] = prod.real();
                    jim[j] = prod.imag();
                }
                
                // integrate the real and imaginary part as independent cubic splines
                gsl_interp_init(spline, inp_.rgrid.data(), jre.data(), rad_.rgrid.size());
                Jre[i] = gsl_interp_eval_integ(spline, inp_.rgrid.data(), jre.data(), 0, rad_.rgrid.back(), acc);
                gsl_interp_init(spline, inp_.rgrid.data(), jim.data(), rad_.rgrid.size());
                Jim[i] = gsl_interp_eval_integ(spline, inp_.rgrid.data(), jim.data(), 0, rad_.rgrid.back(), acc);
            }
            
            // value of the projectile function
            double jf_eval = special::ric_j(ell, kf * r_eval);
            double djf_eval = kf * special::dric_j(ell, kf * r_eval);
            
            // replace the sampled overlaps by a cubic spline and calculate the wronskian with the projectile radial function
            gsl_interp_init(spline, inp_.rgrid.data(), Jre.data(), rad_.rgrid.size());
            double Jre_eval = gsl_interp_eval(spline, inp_.rgrid.data(), Jre.data(), r_eval, acc);
            double dJre_eval = gsl_interp_eval_deriv(spline, inp_.rgrid.data(), Jre.data(), r_eval, acc);
            double Wre = jf_eval * dJre_eval - djf_eval * Jre_eval;
            
            gsl_interp_init(spline, inp_.rgrid.data(), Jim.data(), rad_.rgrid.size());
            double Jim_eval = gsl_interp_eval(spline, inp_.rgrid.data(), Jim.data(), r_eval, acc);
            double dJim_eval = gsl_interp_eval_deriv(spline, inp_.rgrid.data(), Jim.data(), r_eval, acc);
            double Wim = jf_eval * dJim_eval - djf_eval * Jim_eval;
            
            // add this partial T-matrix
            T_matrices[iEpert][istate][fstate][ell] += Complex(Wre,Wim) * 4. * special::constant::pi / kf * std::pow(Complex(0.,1.), -ell)
                * special::ClebschGordan(lf, mf, ell, mi - mf, inp_.L, mi) * special::constant::sqrt_half;
        }
    }
    
    gsl_interp_free(spline);
    gsl_interp_accel_free(acc);
}

void Amplitudes::write ()
{
    for (std::size_t iEpert = 0; iEpert < cmd_.Epert.size(); iEpert++)
    {
        //
        // Write T-matrices to SQL batch file (for use in hex-db).
        //
        
        std::ofstream Tmfile (format("tmat-E%g-%d-%d-%d.sql", inp_.Etot + cmd_.Epert[iEpert], inp_.L, inp_.S, inp_.Pi));
        
        // set exponential format for floating point output
        Tmfile.setf(std::ios_base::scientific);
        
        // write header
        Tmfile << logo("--");
        Tmfile << "-- File generated on " << current_time();
        Tmfile << "--" << std::endl;
        Tmfile << "-- Partial T-matrices for use in the database interface program \"hex-db\"." << std::endl;
        Tmfile << "-- Use for example:" << std::endl;
        Tmfile << "--    > hex-db --new --database hex.db --import <sqlfile> --update" << std::endl;
        Tmfile << "--" << std::endl;
        Tmfile << "BEGIN TRANSACTION;" << std::endl;
        
        for (std::size_t istate = 0; istate < inp_.istates.size(); istate++)
        for (std::size_t fstate = 0; fstate < inp_.fstates.size(); fstate++)
        for (std::size_t ell = 0; ell <= ang_.maxell(); ell++)
        {
            int ni = inp_.istates[istate].n, li = inp_.istates[istate].l, mi = inp_.istates[istate].m;
            int nf = inp_.fstates[fstate].n, lf = inp_.fstates[fstate].l, mf = inp_.fstates[fstate].m;
            
            Complex T = T_matrices[iEpert][istate][fstate][ell];
            
            if (Complex_finite(T) and T != 0.)
            {
                Tmfile  << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                        << ni << ',' << li << ',' << mi << ','
                        << nf << ',' << lf << ',' << mf << ','
                        << inp_.L  << ',' << 0 << ','
                        << inp_.Etot + cmd_.Epert[iEpert] + 1. / (ni * ni) << ',' << ell << ',' 
                        << T.real() << ',' << T.imag() << ");" << std::endl;
            }
        }
        
        Tmfile << "COMMIT;" << std::endl;
    
        //
        // Write cross sections to text file.
        //
        
        std::ofstream csfile (format("cs-E%g.txt", inp_.Etot + cmd_.Epert[iEpert]));
        csfile << logo("#");
        csfile << "# File generated on " << current_time() << "#" << std::endl;
        csfile << "# " << (inp_.S == 0 ? "Singlet" : "Triplet") << " partial cross sections." << std::endl;
        csfile << "#" << std::endl;
        csfile << "# istate / fstate / cross section" << std::endl;
        
        for (std::size_t istate = 0; istate < inp_.istates.size(); istate++)
        for (std::size_t fstate = 0; fstate < inp_.fstates.size(); fstate++)
        {
            int ni = inp_.istates[istate].n, li = inp_.istates[istate].l, mi = inp_.istates[istate].m;
            int nf = inp_.fstates[fstate].n, lf = inp_.fstates[fstate].l, mf = inp_.fstates[fstate].m;
            
            double ki = std::sqrt(inp_.Etot + cmd_.Epert[iEpert] + 1./(ni*ni));
            double kf = std::sqrt(inp_.Etot + cmd_.Epert[iEpert] + 1./(nf*nf));
            
            double cs = T_matrices[iEpert][istate][fstate].sqrnorm() * (2 * inp_.S + 1) / special::pow_int(4. * special::constant::pi, 2) * ki / kf;
            
            csfile << Hydrogen::stateName(ni,li,mi) << '\t' << Hydrogen::stateName(nf,lf,mf) << '\t' << cs << std::endl;
        }
    }
}
