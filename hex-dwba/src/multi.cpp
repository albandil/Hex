/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
 * \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <cstdlib>
#include <sstream>

#include "chebyshev.h"
#include "clenshawcurtis.h"
#include "hydrogen.h"
#include "multi.h"

PhiFunctionDir::PhiFunctionDir (
    HydrogenFunction const & psin, 
    int lam, 
    DistortingPotential const & U, 
    HydrogenFunction const & psi,
    int cblimit
) : psin(psin), psi(psi), Lam(lam), U(U), Diag(psin == psi), 
    Zero(lam == 0 and Diag and (DistortingPotential(psi) == U)),
    Wavenum(psin.getK()), UseFront(false)
{
    if (Zero)
        return;
    
    // compute turning point
    r0 = min(psin.getTurningPoint(), psi.getTurningPoint());
    
    // we need some space around the origin even if there is no turning point
    // FIXME some better way to handle bound-bound integrals instead of the fixed "1."
    if (r0 == 0.)
        r0 = 1.;
    
    // prepare filenames
    std::ostringstream oss;
    oss << "phidir-" << lam << "-"
        << psi.getK()  << "-" << psi.getN()  << "-" << psi.getL()  << "-"
        << psin.getK() << "-" << psin.getN() << "-" << psin.getL() << "~";
    std::string name0 = oss.str() + "f.chb";
    std::string name1 = oss.str() + "1a.chb";
    std::string name2 = oss.str() + "2a.chb";
    
    rArray a, b, a0, b0;
    if (a.hdfload(name1.c_str()) and b.hdfload(name2.c_str()) and (r0 == 0. or a0.hdfload(name0.c_str())))
    {
        Cheb_L = Chebyshev<double,double>(a, -1, 1);
        Cheb_mLm1 = Chebyshev<double,double>(b, -1, 1);
        
        Cheb_L_tail = Cheb_L.tail(1e-10);
        Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10);
        
        if (r0 > 0.)
        {
            Cheb0_L = Chebyshev<double,double>(a0, -1, 1);
            Cheb0_L_tail = Cheb0_L.tail(1e-10);
        }
    }
    else
    {
        // auxiliary variables indicating success when approximating the functions
        bool Cheb_L_conv = false, Cheb_mLm1_conv = false, Cheb0_L_conv = false;
        
        // find expansion of the function in the classically forbidden region
        if (r0 == 0.)
        {
            // no need to approximate anything
            Cheb0_L_conv = true;
        }
        else
        {
            // integral evaluator
            auto inte0 = [&](double x) -> double
            {
                double sharpfactor = 1 / pow(x, 2*Lam+1);
                
                auto integ0 = [&](double r) -> double
                {
                    double ypsi, ypsin; int kpsi, kpsin;
                    std::tie(ypsi,kpsi) = psi.getZeroAsymptotic(r);
                    std::tie(ypsin,kpsin) = psin.getZeroAsymptotic(r);
                    
                    return pow(r, kpsi + kpsin + Lam) * sharpfactor * ypsi * ypsin;
                };
                
                ClenshawCurtis<decltype(integ0), double> cc(integ0);
                return cc.integrate(0., x);
            };
            
            // construct Chebyshev expansion
            for (int N = 4; (cblimit < 0 or N <= cblimit) and not Cheb0_L_conv; N *= 2)
            {
                Cheb0_L.generate(inte0, N, 0, r0);
                
                if ((Cheb0_L_tail = Cheb0_L.tail(1e-10)) < N)
                    Cheb0_L_conv = true;
            }
            
            // check success
            if (not Cheb0_L_conv)
                throw exception ("[PhiFunctionDir] Unable to approximate function in the classically forbidden region.");
            
            Cheb0_L.coeffs().hdfsave(name0.c_str());
        }
        
        // try to find the expansions in the classically allowed region by real integration
        if (psin.getK() < 1.)
        {
            // for low energies this is the only way, so we disable the limit
            tryFullRealChebyshev(cblimit, Cheb_L_conv, Cheb_mLm1_conv);
            
            // check success
            if (Cheb_L_conv and Cheb_mLm1_conv)
            {
                // store precomputed expansions to disk
                Cheb_L.coeffs().hdfsave(name1.c_str());
                Cheb_mLm1.coeffs().hdfsave(name2.c_str());
            }
            else
            {
                throw exception ("[PhiFunctionDir] cblimit (%d) too low for evaluation of the function in the classically allowed region.", cblimit);
            }
        }
        else
        {
            // for high energies we may back away after finite cblimit
            tryFullRealChebyshev(cblimit, Cheb_L_conv, Cheb_mLm1_conv);
            
            // check success
            if (Cheb_L_conv and Cheb_mLm1_conv)
            {
                // store precomputed expansions to disk
                Cheb_L.coeffs().hdfsave(name1.c_str());
                Cheb_mLm1.coeffs().hdfsave(name2.c_str());
                return;
            }
            
            // will approximate the function only on the beginning by the
            // direct Chebyshev expansion, and the rest by approximating
            // the modulation and phase dependences
            tryFrontRealChebyshev(cblimit, Cheb_L_conv, Cheb_mLm1_conv);
            UseFront = true;
        }
    }
}

void PhiFunctionDir::tryFullRealChebyshev (
    int cblimit,
    bool & Cheb_L_conv,
    bool & Cheb_mLm1_conv
){
    std::cout << "try real... " << std::flush;
    
    //
    // integrals in the classically allowed region
    //
    
    auto inte1 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,Lam); };
    auto inte2 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,-Lam-1); };
    
    // integrand compactifications in the classically allowed region
    CompactIntegrand<decltype(inte1),double> compact1(inte1, 0, Inf);
    CompactIntegrand<decltype(inte2),double> compact2(inte2, 0, Inf);
    
    // run the evaluation/convergence loop
    for (int N = 16; (cblimit < 0 or N <= cblimit) and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
    {
        if (not Cheb_L_conv)
        {
            Cheb_L.generate(compact1, N, -1, 1);
            if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
            {
                Cheb_L_conv = true;
                Cheb_L = Cheb_L.integrate(Cheb_L.Integ_Low);
            }
        }
        if (not Cheb_mLm1_conv)
        {
            Cheb_mLm1.generate(compact2, N, -1, 1);
            if ((Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10)) < N)
            {
                Cheb_mLm1_conv = true;
                Cheb_mLm1 = Cheb_mLm1.integrate(Cheb_mLm1.Integ_High);
            }
        }
    }
}

void PhiFunctionDir::tryFrontRealChebyshev (
    int cblimit,
    bool & Cheb_L_conv,
    bool & Cheb_mLm1_conv
){
    std::cout << "try front... " << std::flush;
    
    // integrands
    auto inte1 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,Lam); };
    auto inte2 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,-Lam-1); };
    
    // end of front (ten oscillation periods)
    double end = 20.*M_PI/psin.getK();
    
    // run the evaluation/convergence loop
    for (int N = 16; (cblimit < 0 or N <= cblimit) and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
    {
        if (not Cheb_L_conv)
        {
            Cheb_L.generate(inte1, N, 0., end);
            if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
            {
                Cheb_L_conv = true;
                Cheb_L = Cheb_L.integrate(Cheb_L.Integ_Low);
            }
        }
        if (not Cheb_mLm1_conv)
        {
            Cheb_mLm1.generate(inte2, N, 0., end);
            if ((Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10)) < N)
            {
                Cheb_mLm1_conv = true;
                Cheb_mLm1 = Cheb_mLm1.integrate(Cheb_mLm1.Integ_High);
            }
        }
    }
    
    if (not Cheb_L_conv or not Cheb_mLm1_conv)
        throw exception("[Chebyshev] insufficient cblimit (%d) for front approximation.", cblimit);
    
    // the following code computes A and B from
    //
    //     φ(r₂) = A(r₂)sin(kr₂) + B(r₂)cos(kr₂)
    //
    
    // return amplitude of sine, A(r2)
    auto AmplitudeA = [&](double r2) -> double {
        
        // get left node r_n = (2n+1)*π/(2k)
        int n = floor(psin.getK() * r2 / M_PI - 0.5);
        
        double x2n, x2np1, q1a, q2a, An, Anp1;
        
        // evaluate A at r_n
        x2n = (n + 0.5) * M_PI / psin.getK();
        auto inte1a = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x1/x2n,Lam)/x2n; };
        auto inte2a = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x2n/x1,Lam)/x1; };
        ClenshawCurtis<decltype(inte1a),double> Q1a(inte1a);
        q1a = Q1a.integrate(0.,x2n);
        ClenshawCurtis<decltype(inte2a),double> Q2a(inte2a);
        q2a = Q2a.integrate(x2n,Inf);
        An = (n % 2 == 0 ? q1a + q2a : -q1a - q2a);
        
        // evaluate A at r_n+1
        x2np1 = ((n + 1) + 0.5) * M_PI / psin.getK();
        auto inte1ap = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x1/x2np1,Lam)/x2np1; };
        auto inte2ap = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x2np1/x1,Lam)/x1; };
        ClenshawCurtis<decltype(inte1ap),double> Q1ap(inte1ap);
        q1a = Q1ap.integrate(0.,x2np1);
        ClenshawCurtis<decltype(inte2ap),double> Q2ap(inte2ap);
        q2a = Q2ap.integrate(x2np1,Inf);
        Anp1 = ((n + 1) % 2 == 0 ? q1a + q2a : -q1a - q2a);
        
        // return average
        return (An * (x2np1 - r2) + Anp1 * (r2 - x2n)) / (x2np1 - x2n);
    };
    
    // return amplitude of cosine, B(r2)
    auto AmplitudeB = [&](double r2) -> double {
        
        // get left node r_n = n*π/k
        int n = floor(psin.getK() * r2 / M_PI);
        
        double x2n, x2np1, q1b, q2b, Bn, Bnp1;
        
        // evaluate A at r_n
        x2n = n * M_PI / psin.getK();
        auto inte1b = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x1/x2n,Lam)/x2n; };
        auto inte2b = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x2n/x1,Lam)/x1; };
        ClenshawCurtis<decltype(inte1b),double> Q1b(inte1b);
        q1b = Q1b.integrate(0.,x2n);
        ClenshawCurtis<decltype(inte2b),double> Q2b(inte2b);
        q2b = Q2b.integrate(x2n,Inf);
        Bn = (n % 2 == 0 ? q1b + q2b : -q1b - q2b);
        
        // evaluate A at r_n+1
        x2np1 = (n + 1) * M_PI / psin.getK();
        auto inte1bp = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x1/x2np1,Lam)/x2np1; };
        auto inte2bp = [&](double x1) -> double { return (x1 == 0. or not std::isfinite(x1)) ? 0. : psin(x1)*psi(x1)*pow(x2np1/x1,Lam)/x1; };
        ClenshawCurtis<decltype(inte1bp),double> Q1bp(inte1bp);
        q1b = Q1bp.integrate(0.,x2np1);
        ClenshawCurtis<decltype(inte2bp),double> Q2bp(inte2bp);
        q2b = Q2bp.integrate(x2np1,Inf);
        Bnp1 = ((n + 1) % 2 == 0 ? q1b + q2b : -q1b - q2b);
        
        // return average
        return (Bn * (x2np1 - r2) + Bnp1 * (r2 - x2n)) / (x2np1 - x2n);
    };
    
    // compactifications
    CompactificationR<decltype(AmplitudeA),double> compactA(AmplitudeA);
    CompactificationR<decltype(AmplitudeB),double> compactB(AmplitudeB);
    
    std::cout << "and envelope... ";
    
    // run the evaluation/convergence loop
    bool A_conv = false, B_conv = false;
    for (int N = 16; (cblimit < 0 or N <= cblimit) and not (A_conv and B_conv); N *= 2)
    {
        if (not A_conv)
        {
            ACoeff.generate(compactA, N, compactA.scale(end), 1.);
            if ((ATail = ACoeff.tail(1e-10)) < N)
                A_conv = true;
        }
        if (not B_conv)
        {
            BCoeff.generate(compactB, N, compactB.scale(end), 1.);
            if ((BTail = BCoeff.tail(1e-10)) < N)
                B_conv = true;
        }
        
        std::ostringstream ossA, ossB;
        ossA << "A-" << N << ".dump";
        ossB << "B-" << N << ".dump";
        ACoeff.coeffs().hdfsave(ossA.str().c_str());
        BCoeff.coeffs().hdfsave(ossB.str().c_str());
    }
    if (not A_conv or not B_conv)
    {
        ACoeff.coeffs().hdfsave("ACoeff.dump");
        BCoeff.coeffs().hdfsave("BCoeff.dump");
        throw exception("[Chebyshev] insufficient cblimit (%d) for envelope approximation.", cblimit);
    }
}

void PhiFunctionDir::tryFullComplexChebyshev (
    int cblimit,
    bool & Cheb_L_conv,
    bool & Cheb_mLm1_conv
){
    std::cout << " try complex... " << std::flush;
    
    // expressions to approximate
    auto expr1 = [&](double r2) -> double {
        
        // complex integrand
        auto inte1 = [&](double xi) -> double {
            Complex r1(r2,xi);
            Complex eval = psin.Hplus(r1) * psi(r1) * pow(r1/r2,Lam)/r2;
            return eval.real();
        };
        
        // compute the integral
        return -ClenshawCurtis<decltype(inte1),double>(inte1).integrate(0.,Inf) * pow(r2,Lam+1);
        
    };
    
    auto expr2 = [&](double r2) -> double {
        
        // complex integrand
        auto inte2 = [&](double xi) -> double {
            Complex r1(r2,xi);
            Complex eval = psin.Hplus(r1) * psi(r1) * pow(r2/r1,Lam)/r1;
            return eval.real();
        };
        
        // compute the integral
        return ClenshawCurtis<decltype(inte2),double>(inte2).integrate(0,Inf) * pow(r2,-Lam);
        
    };
    
    // compactifications
    CompactificationR<decltype(expr1),double> compactA(expr1);
    CompactificationR<decltype(expr2),double> compactB(expr2);
    
    // Coulomb wave prefactor which is not included in psin.Hplus()
    double prefactor = sqrt(M_2_PI)/psin.getK();
    
    // run the evaluation/convergence loop
    for (int N = 16; (cblimit < 0 or N <= cblimit) and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
    {
        std::cout << N << " ";
        
        if (not Cheb_L_conv)
        {
            Cheb_L.generate(compactA, N, -1, 1);
            
            if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
            {
                // multiply all coefficients by missing coefficient sqrt(2/π)/k
                Cheb_L = Chebyshev<double,double>(prefactor * Cheb_L.coeffs(), -1, 1);
                Cheb_L_conv = true;
            }
        }
        if (not Cheb_mLm1_conv)
        {
            Cheb_mLm1.generate(compactB, N, -1, 1);
            
            if ((Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10)) < N)
            {
                // multiply all coefficients by missing coefficient sqrt(2/π)/k
                Cheb_mLm1 = Chebyshev<double,double>(prefactor * Cheb_mLm1.coeffs(), -1, 1);
                Cheb_mLm1_conv = true;
            }
        }
    }
}

double PhiFunctionDir::operator() (double r2) const
{
    double t2 = (r2 - 1) / (r2 + 1);
    double tp = 2 * r2 / r0 - 1;
    
    if (Zero)
        return 0.;
    
    if (r2 == 0. and Lam != 0)
        return 0;
    
    else if (r2 == 0 and Lam == 0)
        return Cheb_mLm1.clenshaw(t2, Cheb_mLm1_tail) + Cheb0_L.clenshaw(tp, Cheb0_L_tail);
        
    else if (UseFront)
    {
        if (psin.getK() * r2 < 20. * M_PI)
        {
            // return front Chebyshev-approximated form
            return pow(r2,Lam) * Cheb_mLm1.clenshaw(t2, Cheb_mLm1_tail)
                   + pow(r2, -Lam-1) * Cheb_L.clenshaw(t2, Cheb_L_tail);
        }
        else
        {
            // return envelope approximation
            return ACoeff.clenshaw(t2, ATail) * sin(psin.getK() * r2) +
                   BCoeff.clenshaw(t2, BTail) * cos(psin.getK() * r2);
        }
    }
    else if (r2 >= r0)
    {
        // return full Chebyshev-approximated form in the classically allowed region
        return pow(r2,Lam) * Cheb_mLm1.clenshaw(t2, Cheb_mLm1_tail)
               + pow(r2, -Lam-1) * Cheb_L.clenshaw(t2, Cheb_L_tail);
    }
    else
    {
        // return full Chebyshev-approximated form in the classically forbidden region
        return pow(r2,Lam) * ( Cheb_mLm1.clenshaw(t2, Cheb_mLm1_tail) + Cheb0_L.clenshaw(tp, Cheb0_L_tail) );
    }
}

double PhiFunctionDir::getTurningPoint() const
{
    return r0;
}

std::pair<double,int> PhiFunctionDir::getZeroAsymptotic(double x) const
{
    if (Zero)
        return std::make_pair(0.0, 0);
    
    if (x >= r0)
    {
        // adapt allowed-region solution
        
        double t = (x - 1) / (x + 1);
        return std::make_pair (
            Cheb_mLm1.clenshaw(t, Cheb_mLm1_tail) + (x == 0. ? 0. : pow(x, -2*Lam-1) * Cheb_L.clenshaw(t, Cheb_L_tail)),
            Lam
        );
    }
    else
    {
        // use the prepared forbidden-region solution
        
        double t = (x - 1) / (x + 1);
        double tp = 2 * x / r0 - 1;
        return std::make_pair (
            Cheb_mLm1.clenshaw(t, Cheb_mLm1_tail) + Cheb0_L.clenshaw(tp, Cheb0_L_tail),
            Lam
        );
    }
}
