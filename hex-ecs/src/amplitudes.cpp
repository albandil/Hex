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

#include "arrays.h"
#include "bspline.h"
#include "chebyshev.h"
#include "clenshawcurtis.h"
#include "radial.h"
#include "special.h"
#include "matrix.h"

cArray computeLambda
(
    Bspline const & bspline,
    rArray const & kf, rArray const & ki,
    int maxell, int L, int Spin, int Pi,
    int ni, int li, int mi,
    rArray const & Ei, int lf,
    cArray const & Pf_overlaps,
    std::vector<std::pair<int,int>> const & coupled_states
)
{
    // shorthands
    unsigned Nenergy = kf.size();                // energy count
    Complex const * const t = &(bspline.t(0));   // B-spline knots
    int order   = bspline.order();               // B-spline order
    int Nspline = bspline.Nspline();             // B-spline count
    int Nknot   = bspline.Nknot();               // number of all knots
    int Nreknot = bspline.Nreknot();             // number of real knots
    
    cArray rads(Nenergy * (maxell + 1));
    
    // for all energies, compute the radial factors
    for (unsigned ie = 0; ie < Nenergy; ie++)
    {
        // compose filename of the data file for this solution
        std::ostringstream oss;
        oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
        
        // load the solution
        cArray solution;
        bool solution_exists = false;
        #pragma omp critical
        solution_exists = solution.hdfload(oss.str().c_str());
        if (not solution_exists)
        {
            std::cout << "File \"" << oss.str() << "\" not found." << std::endl;
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
            cArray j_R0 = special::ric_jv(maxell, kf[ie] * eval_r);
            cArray dj_R0 = special::dric_jv(maxell, kf[ie] * eval_r) * kf[ie];
            
            // evaluate B-splines and their derivatives at evaluation radius
            CooMatrix Bspline_R0(Nspline, 1), Dspline_R0(Nspline, 1);
            for (int ispline = 0; ispline < Nspline; ispline++)
            {
                Complex val;
                
                // evaluate B-spline
                val = bspline.bspline(ispline, eval_knot-1, order, eval_r);
                if (val != 0.)
                    Bspline_R0.add(ispline, 0, val);
                
                // evaluate B-spline derivative
                val = bspline.dspline(ispline, eval_knot-1, order, eval_r);
                if (val != 0.)
                    Dspline_R0.add(ispline, 0, val);
            }
            
            // evaluate Wronskians
            CooMatrix Wj[maxell + 1];
            for (int l = 0; l <= maxell; l++)
                Wj[l] = dj_R0[l] * Bspline_R0 - j_R0[l] * Dspline_R0;
                
            // we need "P_overlaps" to have a 'dot' method
            CooMatrix Sp (Nspline, 1, Pf_overlaps.begin());
            
            // compute radial factor
            #pragma omp parallel for
            for (unsigned ill = 0; ill < coupled_states.size(); ill++)
            {
                // use blocks that result in requested final lf
                int l1 = coupled_states[ill].first;
                int l2 = coupled_states[ill].second;
                if (l1 != lf)
                    continue;
                
                // get correct solution (for this ang. mom.)
                cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
                
                rads[ie * (maxell + 1) + l2] += Sp.transpose().dot(PsiSc).dot(Wj[l2].todense()).todense()[0] / double(samples);
            }
        }
    }
    
    return rads;
}

Chebyshev<double,Complex> fcheb (Bspline const & bspline, cArrayView const & PsiSc, double kmax, int l1, int l2)
{
    // shorthands
    Complex const * const t = &(bspline.t(0));   // B-spline knots
    int Nspline = bspline.Nspline();             // number of real knots
    int Nreknot = bspline.Nreknot();             // number of real knots
    int order   = bspline.order();               // B-spline order
    
    // determine evaluation radius
    char const * HEX_RHO = getenv("HEX_RHO");
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
            
            /// DEBUG
            if (err1 != GSL_SUCCESS or err2 != GSL_SUCCESS)
            {
                std::cerr << "Errors while evaluating Coulomb function:\n";
                std::cerr << "\terr1 = " << err1 << "\n";
                std::cerr << "\terr2 = " << err2 << "\n";
                exit(-1);
            }
            ///
            
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
                std::cerr << "F1F2 = " << F1F2 << "\n";
            if (not std::isfinite(std::abs(ddrho_Psi)))
                std::cerr << "ddrho_Psi = " << ddrho_Psi << "\n";
            if (not std::isfinite(std::abs(Psi)))
                std::cerr << "Psi = " << Psi << "\n";
            if (not std::isfinite(ddrho_F1F2))
                std::cerr << "ddrho_F1F2 = " << ddrho_F1F2 << "\n";
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
            throw exception("ERROR: Non-convergent Chebyshev expansion.");
    }
    
    return CB;
}

cArrays computeXi
(
    Bspline const & bspline, int maxell, int L, int Spin, int Pi, int ni, int li, int mi, 
    rArray const & Ei, rArray & ics, std::vector<std::pair<int,int>> const & coupled_states
)
{
    // resize and clear the output storage for integral cross sections
    ics.resize(Ei.size());
    ics.clear();
    
    // array of Chebyshev expansions for every Ei and angular state
    cArrays results (coupled_states.size());
    
    // B-spline count
    int Nspline = bspline.Nspline();
    
    // for all energies
    for (size_t ie = 0; ie < Ei.size(); ie++)
    {
        // compose filename of the data file for this solution
        std::ostringstream oss;
        oss << "psi-" << L << "-" << Spin << "-" << Pi << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
        
        // load the solution
        cArray solution;
        # pragma omp critical
        if (not solution.hdfload(oss.str().c_str()))
            throw exception ("Can't open the solution file \"%s\"!", oss.str().c_str());
        
        // maximal available momentum
        double kmax = sqrt(Ei[ie] - 1./(ni*ni));
        
        // for all angular states ???: (triangle ℓ₂ ≤ ℓ₁)
        for (unsigned ill = 0; ill < coupled_states.size(); ill++)
        {
            int l1 = coupled_states[ill].first;
            int l2 = coupled_states[ill].second;
            
            std::cout << "\tEi[" << ie << "] = " << Ei[ie] << ", l1 = " << l1 << ", l2 = " << l2 << "\n";
            
            // create subset of the solution
            cArrayView PsiSc (solution, ill * Nspline * Nspline, Nspline * Nspline);
            
            // compute new ionization amplitude
            Chebyshev<double,Complex> CB = fcheb (bspline, PsiSc, kmax, l1, l2);
            results[ill] = CB.coeffs();
            
            // integrate the expansion
            int tail = CB.tail(1e-10);
            int n;
            auto fsqr = [&](double beta) -> double { return sqrabs(CB.clenshaw(kmax*sin(beta), tail)); };
            ClenshawCurtis<decltype(fsqr),double> integrator(fsqr);
            double cs = integrator.integrate(0, special::constant::pi_quart, &n) / sqrt(Ei[ie]);
            
            std::cout << "\t\t- contrib to ics: " << cs << " (" << n << " evaluations)\n";
            ics[ie] += cs;
        }
    }
    
    return results;
}
