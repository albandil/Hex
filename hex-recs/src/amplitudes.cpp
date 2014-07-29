/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <regex>
#include <vector>

#include <fftw3.h>

#include "amplitudes.h"
#include "arrays.h"
#include "bspline.h"
#include "chebyshev.h"
#include "clenshawcurtis.h"
#include "hydrogen.h"
#include "radial.h"
#include "special.h"
#include "matrix.h"

void Amplitudes::extract ()
{
    std::cout << std::endl << "Extracting T-matrices" << std::endl;
    
    // for all initial and final states
    for (auto instate  : inp_.instates)
    for (auto outstate : inp_.outstates)
    {
        // get initial quantum numbers
        int ni = std::get<0>(instate);
        int li = std::get<1>(instate);
        int two_ji = std::get<2>(instate); double ji = 0.5 * two_ji;
        int two_mi = std::get<3>(instate); double mi = 0.5 * two_mi;
        
        // get final quantum numbers
        int nf = std::get<0>(outstate);
        int lf = std::get<1>(outstate);
        int two_jf = std::get<2>(outstate); double jf = 0.5 * two_jf;
        
        // non-relativistic total energies
        rArray ETot = inp_.Ei - 1./(ni*ni) /*- mi * inp.B*/;
        
        // add spin-orbital energy
        if (li > 0)
            ETot += 0.5 * special::constant::alpha_sqr * (ji*(ji+1)-li*(li+1)-0.75) * Hydrogen::moment(-3,ni,li);
        
        if (nf > 0)
        {
            //
            // Discrete transition
            //
            
            std::cout << "\texc: -> nf = " << nf << ", lf = " << lf << ", jf = " << jf << ", mf = *" << std::flush;
            
            // precompute hydrogen function overlaps
            cArray Pf_overlaps = rad_.overlapP(nf,lf,weightEndDamp(bspline_));
            
            // compute radial integrals
            for (int two_mf = -two_jf; two_mf <= two_jf; two_mf++)
            {
                // final atomic energies
                rArray Ef = -1./(nf*nf) /*- mf * inp.B*/;
                
                // add spin-orbital energy
                if (lf > 0)
                    Ef += 0.5 * special::constant::alpha_sqr * (jf*(jf+1)-lf*(lf+1)-0.75) * Hydrogen::moment(-3,nf,lf);
                
                // final projectile momenta
                rArray kf = sqrt(ETot - Ef);
                
                // setup transition name
                auto transition = std::make_tuple(ni,li,two_ji,two_mi,nf,lf,two_jf,two_mf);
                
                // compute Λ for this transition
                computeLambda_(transition, kf);
                
                // compute T-matrices for this transition
                computeTmat_(transition, kf);
                
                // compute cross sections for this transition
                computeSigma_(transition, kf);
            }
            
            std::cout << " ok" << std::endl;
        }
        else
        {
            //
            // Ionization
            //
            
            throw exception ("Ionization not implemented yet.");
        }
    }
}

std::map<std::tuple<int,int,int>,cArray> Amplitudes::computeLambda_
(
    std::tuple<int,int,int,int,int,int,int,int> transition,
    rArray const & kf
)
{
    // shorthands
    unsigned Nenergy = kf.size();                // energy count
    Complex const * const t = &(bspline_.t(0));  // B-spline knots
    int order   = bspline_.order();              // B-spline order
    int Nspline = bspline_.Nspline();            // B-spline count
    int Nknot   = bspline_.Nknot();              // number of all knots
    int Nreknot = bspline_.Nreknot();            // number of real knots
    
    // get transition information
    int ni = std::get<0>(transition);      int nf = std::get<4>(transition);
    int li = std::get<1>(transition);      int lf = std::get<5>(transition);
    int two_ji = std::get<2>(transition);  int two_jf = std::get<6>(transition);
    int two_mi = std::get<3>(transition);  int two_mf = std::get<6>(transition);
    
    // compute final hydrogen orbital overlaps with B-spline basis
    cArray Pf_overlaps = rad_.overlapP(nf, lf, weightEndDamp(bspline_));
    
    // output array
    std::map<std::tuple<int,int,int>,cArray> Lambda;
    
    // for all energies, compute the radial factors
    for (unsigned ie = 0; ie < Nenergy; ie++)
    {
        // compose filename of the data file for this solution
        SolutionIO reader (inp_.J, ni, li, two_ji, two_mi, inp_.Ei[ie]);
        
        // load the solution
        cArray solution;
        if (not reader.load(solution))
        {
            std::cout << "File \"" << reader.name() << "\" not found." << std::endl;
            continue;
        }
        
        // The cross section oscillates, so we will do some averaging
        // As recommended by Bartlett, we will compute several amplitudes
        // separated by π/(n*kf[ie]) near the R₀ turning point.
        double wavelength = special::constant::pi / kf[ie];
        char const * HEX_RHO = getenv("HEX_RHO");
        char const * HEX_SAMPLES = getenv("HEX_SAMPLES");
        int samples = (HEX_SAMPLES == nullptr) ? 10 : atoi(HEX_SAMPLES);
        double R0 = (HEX_RHO == nullptr) ? t[Nreknot - 1].real() : atof(HEX_RHO);
        
        // skip impact energies with undefined outgoing momentum
        if (std::isnan(kf[ie]))
            continue;
        
        for (int n = 1; n <= samples; n++)
        {
            // this is the evaluation point
            double eval_r = R0 - wavelength * n / samples;
            
            // determine knot
            int eval_knot = std::lower_bound
            (
                t,
                t + Nknot,
                Complex(eval_r, 0.),
                [](Complex a, Complex b) -> bool { return a.real() < b.real(); }
            ) - t;
            
            // evaluate j and dj at far radius for all angular momenta up to maxell
            cArray j_R0 = special::ric_jv(inp_.maxell, kf[ie] * eval_r);
            cArray dj_R0 = special::dric_jv(inp_.maxell, kf[ie] * eval_r) * kf[ie];
            
            // evaluate B-splines and their derivatives at evaluation radius
            CooMatrix Bspline_R0(Nspline, 1), Dspline_R0(Nspline, 1);
            for (int ispline = 0; ispline < Nspline; ispline++)
            {
                Complex val;
                
                // evaluate B-spline
                val = bspline_.bspline(ispline, eval_knot-1, order, eval_r);
                if (val != 0.)
                    Bspline_R0.add(ispline, 0, val);
                
                // evaluate B-spline derivative
                val = bspline_.dspline(ispline, eval_knot-1, order, eval_r);
                if (val != 0.)
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
                if (ang_[ill].l1 != lf)
                    continue;
                
                // get angular momentum
                int ell = ang_[ill].l2;
                
                // get correct solution (for this ang. mom.)
                cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
                
                // initialize the storage (if needed)
                auto index = std::make_tuple(ang_[ill].L,ang_[ill].S,ell);
                if (Lambda.find(index) == Lambda.end())
                    Lambda[index] = cArray(Nenergy);
                
                // update value in the storage
                Lambda[index][ie] += Sp.transpose().dot(PsiSc).dot(Wj[ell].todense()).todense()[0] / double(samples);
            }
        }
    }
    
    return Lambda;
}

std::map<std::tuple<int,int>,cArray> Amplitudes::computeTmat_
(
    std::tuple<int,int,int,int,int,int,int,int> transition,
    rArray const & kf
)
{
    // output array
    std::map<std::tuple<int,int>,cArray> Tmat;
    
    // get transition information
    int ni = std::get<0>(transition);      int nf = std::get<4>(transition);
    int li = std::get<1>(transition);      int lf = std::get<5>(transition);
    int two_ji = std::get<2>(transition);  int two_jf = std::get<6>(transition);
    int two_mi = std::get<3>(transition);  int two_mf = std::get<6>(transition);
    
    // for all j'
    for (int two_jp = std::abs(2*inp_.J-two_jf); two_jp <= 2*inp_.J+two_jf; two_jp += 2)
    {
        double C1 = gsl_sf_coupling_3j
        (
            two_jf,  two_jp,            2*inp_.J,
            two_mf,  inp_.two_M-two_mf, -inp_.two_M
        ) * std::sqrt(2*inp_.J+1);
        
        // for all l'
        for (int lp = std::abs(two_jp-1)/2; lp <= (two_jp+1)/2; lp++)
        {
            double C2 = gsl_sf_coupling_3j
            (
                2*lp,                      1,       two_jp,
                inp_.two_M-two_mf-two_sf,  two_sf,  inp_.two_M-two_mf
            ) * std::sqrt(two_jp+1);
            
            // sum over L and S
            for (auto const & Lambda : Lambda_LSlp[transition])
            {
                int L = std::get<0>(Lambda.first);
                int S = std::get<1>(Lambda.first);
                int ll = std::get<2>(Lambda.first);
                
                // skip terms that do not contribute to (ll = ) lp
                if (ll != lp)
                    continue;
                
                // compute 9j coef
                double W = gsl_sf_coupling_9j
                (
                    2*lf, 1, two_jf,
                    2*lp, 1, two_jp,
                    2*L, 2*S, 2*inp_.J
                );
                
                // initialize the storage (if needed)
                auto index = std::make_tuple(lf,two_jf);
                if (Tmat.find(index) == Tmat.end())
                    Tmat[index] = cArray(inp_.Ei.size());
                
                // add contribution to the T-matrix
                Tmat[index] += C1 * C2 * std::sqrt((two_jf+1)*(two_jp+1)*(2*S+1)*(2*L+1)) * W * Lambda.second / kf;
            }
        }
    }
    
    return Tmat;
}

rArray Amplitudes::computeSigma_
(
    std::tuple<int,int,int,int,int,int,int,int> transition,
    rArray const & kf
)
{
    // get transition information
    int ni = std::get<0>(transition);      int nf = std::get<4>(transition);
    int li = std::get<1>(transition);      int lf = std::get<5>(transition);
    int two_ji = std::get<2>(transition);  int two_jf = std::get<6>(transition);
    int two_mi = std::get<3>(transition);  int two_mf = std::get<6>(transition);
    
    // output array
    rArray sigma (inp_.Ei.size());
    
    // for all summation indices
    for (int two_jp = std::abs(2*inp_.J-two_jf); two_jp <= 2*inp_.J+two_jf; two_jp += 2)
    for (int two_jpp= std::abs(2*inp_.J-two_jf); two_jpp<= 2*inp_.J+two_jf; two_jpp+= 2)
    for (int lp = std::abs(two_jp-1)/2; lp <= (two_jp+1)/2; lp ++)
    for (int lpp= std::abs(two_jp-1)/2; lpp<= (two_jp+1)/2; lpp++)
    if (lp == lpp)
    {
        rArray Re_f_ell = -realpart(Tmat_jplp[transition][std::make_tuple(lp,two_jp)]) / special::constant::two_pi;
        rArray Im_f_ell = -imagpart(Tmat_jplp[transition][std::make_tuple(lp,two_jpp)]) / special::constant::two_pi;
        sigma += kf / inp_.ki * (Re_f_ell * Re_f_ell + Im_f_ell * Im_f_ell);
    }
    
    return sigma;
}

void Amplitudes::writeSQL_files ()
{
    // open file
    std::ofstream fsql (format("tmat-J%d-M%d.sql", inp_.J, inp_.two_M));
    fsql << "BEGIN TRANSACTION;" << std::endl;
    
    // for all transitions
    for (auto data : Tmat_jplp)
    {
        // get transition information
        int ni = std::get<0>(data.first);      int nf = std::get<4>(data.first);
        int li = std::get<1>(data.first);      int lf = std::get<5>(data.first);
        int two_ji = std::get<2>(data.first);  int two_jf = std::get<6>(data.first);
        int two_mi = std::get<3>(data.first);  int two_mf = std::get<6>(data.first);
        
        // for all expansion terms
        for (auto t : data.second)
        {
            // get expansion indices
            int lp = std::get<0>(t.first);
            int two_jp = std::get<1>(t.first);
            
            // get energy-sampled data array
            cArray const & T = t.second;
            
            // for all energies
            for (unsigned i = 0; i < T.size(); i++)
            {
                if (Complex_finite(T[i]) and T[i] != 0.)
                {
                    fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
                         << ni << "," << li << "," << two_ji << "," << two_mi << ","
                         << nf << "," << lf << "," << two_jf << "," << two_mf << ","
                         << inp_.J  << "," << inp_.two_M << ","
                         << inp_.Ei[i] << "," << lp << "," << two_jp << "," 
                         << T[i].real() << "," << T[i].imag() << ","
                         << "0, 0);" << std::endl;
                }
            }
        }
    }
    
    // close file
    fsql << "COMMIT;" << std::endl;
    fsql.close();
}

void Amplitudes::writeICS_files ()
{
    // open file
    std::ofstream fout (format("ics-J%d-M%d.dat", inp_.J, inp_.two_M));
    
    // print table header
    fout << "#E[Ry]\t";
    for (auto data : sigma)
    {
        fout << format
        (
            "%s-%s\t",
            Hydrogen::stateRName
            (
                std::get<0>(data.first),std::get<1>(data.first),std::get<2>(data.first),std::get<3>(data.first)
            ).c_str(),
            Hydrogen::stateRName
            (
                std::get<4>(data.first),std::get<5>(data.first),std::get<6>(data.first),std::get<7>(data.first)
            ).c_str()
        );
    }
    fout << std::endl;
    
    // print data (cross sections)
    for (unsigned ie = 0; ie < inp_.Ei.size(); ie++)
    {
        fout << inp_.Ei[ie] << '\t';
        for (auto data : sigma)
        {
            if (std::isfinite(data.second[ie]))
                fout << data.second[ie] << '\t';
            else
                fout << 0.0 << '\t';
        }
        fout << std::endl;
    }
    
    // close file
    fout.close();
}
