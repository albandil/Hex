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

void extract
(
    CommandLine const & cmd,
    InputFile const & inp,
    AngularBasis const & ang,
    RadialBasis const & rad,
    double Epert,
    Complex const * psi
)
{
    gsl_interp * spline = gsl_interp_alloc(gsl_interp_cspline, rad.rgrid.size());
    gsl_interp_accel * acc = gsl_interp_accel_alloc();
    
    cArrays T_matrices (inp.istates.size() * inp.fstates.size());
    rArray cross_sections (inp.istates.size() * inp.fstates.size());
    
    for (std::size_t istate = 0; istate < inp.istates.size(); istate++)
    {
        // initial state quantum numbers
        int ni = inp.fstates[istate].n;
        //int li = inp.fstates[istate].l; // - not used
        int mi = inp.fstates[istate].m;
        
        // initial projectile  momentum
        double ki = std::sqrt(inp.Etot + Epert + 1./(ni*ni));
        
        for (std::size_t fstate = 0; fstate < inp.fstates.size(); fstate++)
        {
            // final state euantum numbers
            int nf = inp.fstates[fstate].n;
            int lf = inp.fstates[fstate].l;
            int mf = inp.fstates[fstate].m;
            
            // final projectile momentum
            double kf = std::sqrt(inp.Etot + Epert + 1./(nf*nf));
            
            // partial T-matrices
            T_matrices[istate * inp.fstates.size() + fstate] = cArray(ang.maxell() + 1);
            
            // skip in-accessible states
            if (not std::isfinite(kf))
                continue;
            
            // evaluate the final hydrogen function in grid points
            rArray P (rad.rgrid.size());
            for (std::size_t i = 0; i < rad.rgrid.size(); i++)
                P[i] = Hydrogen::P(nf, lf, rad.rgrid[i], inp.Z);
            
            // for all contributing angular blocks
            for (std::size_t iblock = 0; iblock < ang.size(); iblock++)
            if ((unsigned)inp.fstates[fstate].l == ang.state(iblock).first)
            {
                rArray jre (rad.rgrid.size()), jim (rad.rgrid.size()), Jre (rad.rgrid.size()), Jim (rad.rgrid.size());
                
                // extraction radius
                double r_eval = rad.rgrid.back(1);
                
                // partial wave (projectile final angular momentum)
                unsigned ell = ang.state(iblock).second;
                
                // calculate overlaps of all products and the final radial state
                for (std::size_t i = 0; i < rad.rgrid.size(); i++)
                {
                    // project column to radial function
                    for (std::size_t j = 0; j < rad.rgrid.size(); j++)
                    {
                        Complex prod = P[j] * psi[(iblock * rad.Npts + i) * rad.Npts + j];
                        jre[j] = prod.real();
                        jim[j] = prod.imag();
                    }
                    
                    // integrate the real and imaginary part as independent cubic splines
                    gsl_interp_init(spline, inp.rgrid.data(), jre.data(), rad.rgrid.size());
                    Jre[i] = gsl_interp_eval_integ(spline, inp.rgrid.data(), jre.data(), 0, rad.rgrid.back(), acc);
                    gsl_interp_init(spline, inp.rgrid.data(), jim.data(), rad.rgrid.size());
                    Jim[i] = gsl_interp_eval_integ(spline, inp.rgrid.data(), jim.data(), 0, rad.rgrid.back(), acc);
                }
                
                // value of the projectile function
                double jf_eval = special::ric_j(ell, kf * r_eval);
                double djf_eval = kf * special::dric_j(ell, kf * r_eval);
                
                // replace the sampled overlaps by a cubic spline and calculate the wronskian with the projectile radial function
                gsl_interp_init(spline, inp.rgrid.data(), Jre.data(), rad.rgrid.size());
                double Jre_eval = gsl_interp_eval(spline, inp.rgrid.data(), Jre.data(), r_eval, acc);
                double dJre_eval = gsl_interp_eval_deriv(spline, inp.rgrid.data(), Jre.data(), r_eval, acc);
                double Wre = jf_eval * dJre_eval - djf_eval * Jre_eval;
                
                gsl_interp_init(spline, inp.rgrid.data(), Jim.data(), rad.rgrid.size());
                double Jim_eval = gsl_interp_eval(spline, inp.rgrid.data(), Jim.data(), r_eval, acc);
                double dJim_eval = gsl_interp_eval_deriv(spline, inp.rgrid.data(), Jim.data(), r_eval, acc);
                double Wim = jf_eval * dJim_eval - djf_eval * Jim_eval;
                
                // add this partial T-matrix
                T_matrices[istate * inp.fstates.size() + fstate][ell] += Complex(Wre,Wim) * 4. * special::constant::pi / kf * std::pow(Complex(0.,1.), -ell)
                    * special::ClebschGordan(lf, mf, ell, mi - mf, inp.L, mi) * special::constant::sqrt_half;
            }
            
            // calculate cross sections for this L, S and Pi
            cross_sections[istate * inp.fstates.size() + fstate] = (T_matrices[istate * inp.fstates.size() + fstate] / special::constant::two_pi).sqrnorm() * 0.25 * (2 * inp.S + 1) * kf / ki;
        }
    }
    
    gsl_interp_free(spline);
    gsl_interp_accel_free(acc);
    
    //
    // Write T-matrices to SQL batch file (for use in hex-db).
    //
    
    std::ofstream Tmfile (format("tmat-E%g-%d-%d-%d.sql", inp.Etot + Epert, inp.L, inp.S, inp.Pi));
    
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
    
    for (std::size_t istate = 0; istate < inp.istates.size(); istate++)
    {
        for (std::size_t fstate = 0; fstate < inp.fstates.size(); fstate++)
        {
            for (std::size_t ell = 0; ell <= ang.maxell(); ell++)
            {
                int ni = inp.istates[istate].n, li = inp.istates[istate].l, mi = inp.istates[istate].m;
                int nf = inp.fstates[fstate].n, lf = inp.fstates[fstate].l, mf = inp.fstates[fstate].m;
                
                Complex T = T_matrices[istate * inp.fstates.size() + fstate][ell];
                
                if (Complex_finite(T) and T != 0.)
                {
                    Tmfile  << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                            << ni << "," << li << "," << mi << ","
                            << nf << "," << lf << "," << mf << ","
                            << inp.L  << "," << 0 << ","
                            << inp.Etot + Epert + 1. / (ni * ni) << "," << ell << "," 
                            << T.real() << "," << T.imag() << ");" << std::endl;
                }
            }
        }
    }
    
    Tmfile << "COMMIT;" << std::endl;
    
    //
    // Write cross sections to text file.
    //
    
    std::ofstream csfile (format("cs-%g.txt", inp.Etot + Epert));
    csfile << logo("#");
    csfile << "# File generated on " << current_time() << "#" << std::endl;
    csfile << "# " << (inp.S == 0 ? "Singlet" : "Triplet") << " partial cross sections." << std::endl;
    csfile << "#" << std::endl;
    csfile << "# istate / fstate / cross section" << std::endl;
    
    for (std::size_t istate = 0; istate < inp.istates.size(); istate++)
    {
        for (std::size_t fstate = 0; fstate < inp.fstates.size(); fstate++)
        {
            int ni = inp.istates[istate].n, li = inp.istates[istate].l, mi = inp.istates[istate].m;
            int nf = inp.fstates[fstate].n, lf = inp.fstates[fstate].l, mf = inp.fstates[fstate].m;
            
            double cs = cross_sections[istate * inp.fstates.size() + fstate];
            
            csfile << Hydrogen::stateName(ni,li,mi) << "\t" << Hydrogen::stateName(nf,lf,mf) << "\t" << cs << std::endl;
        }
    }
}
