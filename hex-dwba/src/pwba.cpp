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

#include "arrays.h"
#include "complex.h"
#include "gausskronrod.h"
#include "hydrogen.h"
#include "pwba.h"
#include "special.h"
#include "symbolic.h"

using special::constant::pi;
using special::constant::Inf;

void pwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    bool direct, bool exchange
)
{
    // allocate memory
    Tdir.resize((2*Li+1)*(2*Lf+1));
    Texc.resize((2*Li+1)*(2*Lf+1));
    
    // for all outgoing partial waves
    for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
    {
        // add new T-matrix for this outgoing partial wave
        for (cArray & T : Tdir)
            T.push_back(0.);
        for (cArray & T : Texc)
            T.push_back(0.);
        
        // for all incoming partial waves
        for (int li = std::abs(Li - L); li <= Li + L; li++)
        {
            // conserve parity
            if ((lf + Lf) % 2 != (li + Li) % 2)
                continue;
            
            //
            // compute direct contribution
            //
            
            // for all multipoles
            for (int lam = std::max(std::abs(Lf-Li), std::abs(lf-li)); direct and lam <= std::min(Li+Lf, lf+li); lam++)
            {
                // compute the needed radial integrals
                double Vdir = PWBA1::compute_Idir (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
                
                // compute complex prefactor
                Complex prefactor = std::pow(4*pi,2)/(ki*kf) * std::pow(Complex(0.,1.),li-lf) * std::sqrt((2*li+1)/(4*pi));
                
                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
                    
                    // compute angular integrals (Gaunt coefficients)
                    double ang = special::ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi)
                               * special::ClebschGordan(Li,Mi,li,0,L,Mi)
                               * special::computef(lam,Lf,lf,Li,li,L);
                    
                    // add the T-matrix contributions
                    Tdir[idx][lf-std::abs(Lf - L)] += prefactor * ang * Vdir;
                }
            }
            
            //
            // compute exchange contribution
            //
            
            // for all multipoles
            for (int lam = std::max(std::abs(lf-Li), std::abs(Lf-li)); exchange and lam <= std::min(lf+Li, Lf+li); lam++)
            {
                // compute the needed radial integrals
                double Vexc = PWBA1::compute_Iexc (li, lf, lam, Ni, Li, ki, Nf, Lf, kf);
                
                // compute complex prefactor
                Complex prefactor = std::pow(4*pi,2)/(ki*kf)*std::pow(Complex(0.,1.),li-lf)*std::sqrt((2*li+1)/(4*pi));
                
                // for all projections of the initial/final angular momentum
                for (int Mi = -Li; Mi <= Li; Mi++)
                for (int Mf = -Lf; Mf <= Lf; Mf++)
                {
                    // compute index in the array of T-matrices
                    int idx = (Mi + Li)*(2*Lf + 1) + Mf + Lf;
                    
                    // compute angular integrals (Gaunt coefficients)
                    double ang = special::ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi)
                               * special::ClebschGordan(Li,Mi,li,0,L,Mi)
                               * special::computef(lam,lf,Lf,Li,li,L);
                    
                    // add the T-matrix contributions
                    Texc[idx][lf-std::abs(Lf - L)] += prefactor * ang * Vexc;
                }
            }
        }
    }
}

template <class Functor, class Integrator>
class BesselNodeIntegrator1D
{
    private:
        
        Integrator Q_;
        
        double k_;
        int l_;
        
        double result_;
        bool ok_;
        std::string status_;
        
        double epsabs_;
        double epsrel_;
        int limit_;
        
    
    public:
        
        BesselNodeIntegrator1D (Functor f, double k, int l)
            : Q_(f), k_(k), l_(l), result_(0), ok_(true), status_(),
              epsabs_(1e-8), epsrel_(1e-5), limit_(100)
        {
        }
        
        double result () const { return result_; }
        bool ok () const { return ok_; }
        std::string const & status () const { return status_; }
        
        double epsabs () const { return epsabs_; }
        void setEpsAbs (double eps) { epsabs_ = eps; }
        
        double epsrel () const { return epsrel_; }
        void setEpsRel (double eps) { epsrel_ = eps; }
        
        int limit () const { return limit_; }
        void setLimit (double n) { limit_ = n; }
        
        bool integrate (double a, double b)
        {
            // set parcel integrator accuracy to one order more
            Q_.setEpsAbs(0.1 * epsabs_);
            Q_.setEpsRel(0.1 * epsrel_);
            
            // overall integral
            double integral = 0;
            
            // end of previous integration parcel
            double prevR = 0;
            
            // for all integration parcels (nodes of the Bessel function)
            for (int inode = 1; inode < limit_; inode++)
            {
                // get integration bounds
                gsl_sf_result res;
                int err = gsl_sf_bessel_zero_Jnu_e(l_+0.5,inode,&res);
                if (err != GSL_SUCCESS)
                {
                    throw exception
                    (
                        "Cannot find %d-th root of the Bessel function j[%d](%g r) -- %s.",
                        inode, l_, k_, gsl_strerror(err)
                    );
                }
                double rmin = prevR;
                double rmax = res.val / k_;
                
                // skip intervals that are below the lower limit
                if (rmax < a)
                    continue;
                
                // skip intervals that are above the upper limit
                if (rmin > b)
                    break;
                
                // shrink integration interval, if necessary
                rmin = std::max (a, rmin);
                rmax = std::min (rmax, b);
                
                // integrate and check success
                if (not Q_.integrate(rmin,rmax))
                {
                    ok_ = false;
                    status_ = Q_.status();
                    result_ = integral;
                    return ok_;
                }
                
                // update result
                integral += Q_.result();
                prevR = rmax;
                
                // check convergence
                if (std::abs(Q_.result()) < epsabs_ or std::abs(Q_.result()) < epsrel_ * std::abs(integral))
                    break;
            }
            
            ok_ = true;
            result_ = integral;
            status_ = "";
            return ok_;
        }
};

double PWBA1::compute_Idir (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    std::cout << format
    (
        "Precompute Idir\n"
        "\tlambda = %d\n"
        "\tNi = %d, Li = %d, ki = %g, li = %d\n"
        "\tNf = %d, Lf = %d, kf = %g, lf = %d\n",
        lambda, Ni, Li, ki, li, Nf, Lf, kf, lf
    );
    
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf](double r2) -> double
        {
            // inner integrand
            auto iintegrand = [Ni,Li,Nf,Lf,r2](double r1) -> double { return Hydrogen::P(Ni,Li,r1) * Hydrogen::P(Nf,Lf,r1) * (1./r1 - 1./r2); };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            Qi.setEpsAbs(0);
            
            // integrate and check success
            if (not Qi.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Idir (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }
            
            return Qi.result() * special::ric_j(li,ki*r2) * special::ric_j(lf,kf*r2);
        };
        
        // which Bessel function oscillates slowlier ?
        int    l = (ki < kf ? li : lf);
        double k = (ki < kf ? ki : kf);
        
        // outer integrator
        BesselNodeIntegrator1D<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.setEpsAbs(0);
        R.integrate(0,Inf);
        
        std::cout << "\tIdir = " << R.result() << std::endl;
        return R.result();
    }
    else
    {
        // compute normalization factors of the generalized Laguerre polynomials
        double Normi = std::sqrt(std::pow(2./Ni,3) * gsl_sf_fact(Ni-Li-1) / (2. * Ni * gsl_sf_fact(Ni+Li)));
        double Normf = std::sqrt(std::pow(2./Nf,3) * gsl_sf_fact(Nf-Lf-1) / (2. * Nf * gsl_sf_fact(Nf+Lf)));
        
        // compute the polynomials
        symbolic::poly Lagi = symbolic::GeneralizedLaguerre (Ni-Li-1, 2*Li+1);
        symbolic::poly Lagf = symbolic::GeneralizedLaguerre (Nf-Lf-1, 2*Lf+1);
        
        // multiply by the angular factor
        for (symbolic::term & pi : Lagi) pi.a += Li + 1; Normi *= std::pow(2./Ni,Li);
        for (symbolic::term & pf : Lagf) pf.a += Lf + 1; Normf *= std::pow(2./Nf,Lf);
        
        // compute the product of the polynomials
        symbolic::poly PP = Lagi * Lagf;
        
        // factor in the argument of the exponential
        double c = 1./Ni + 1./Nf;
        
        // outer integrand evaluated at "r"
        auto integrand = [PP,c,lambda,Normi,Normf,ki,kf,li,lf](double r) -> double
        {
            // inner integral
            double integral = 0;
            
            // integrate term by term
            for (symbolic::term const & p : PP)
            {
                // compute the high integral
                gsl_sf_result res;
                int err_high = gsl_sf_gamma_inc_e (p.a - lambda, c * r, &res);
                if (err_high != GSL_SUCCESS and err_high != GSL_EUNDRFLW)
                {
                    throw exception
                    (
                        "Unable to evaluate incomplete gamma-function Gamma(%d,%g) - %s.",
                        p.a - lambda, c * r, gsl_strerror(err_high)
                    );
                }
                double int_high = 0;
                if (err_high != GSL_EUNDRFLW)
                    int_high = gsl_sf_pow_int(c*r,lambda) * res.val;
                
                // compute the low integral
                int err_low = gsl_sf_gamma_inc_P_e (p.a + lambda + 1, c * r, &res);
                double scale = gsl_sf_gamma (p.a + lambda + 1);
                if (err_low != GSL_SUCCESS)
                {
                    throw exception
                    (
                        "Unable to evaluate scaled complementary incomplete gamma-function P(%d,%g) - %s.",
                        p.a + lambda + 1, c * r, gsl_strerror(err_low)
                    );
                }
                double int_low = gsl_sf_pow_int(c*r,-lambda-1) * res.val * scale;
                
                // sum both contributions
                integral += (int_low + int_high) * symbolic::double_approx(p.kr) / gsl_sf_pow_int(c,p.a);
            }
            
            return Normi * Normf * integral * special::ric_j(li,ki*r) * special::ric_j(lf,kf*r);
        };
        
        // which Bessel function oscillates slowlier ?
        int    l = std::max(li, lf);
        double k = (ki == kf ? ki : std::abs(ki - kf));
        
        // outer integrator
        BesselNodeIntegrator1D<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.integrate (0,Inf);
        std::cout << "\tIdir = " << R.result() << std::endl;
        
        return R.result();
    }
}

double PWBA1::compute_Iexc (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
{
    std::cout << format
    (
        "Precompute Iexc\n"
        "\tlambda = %d\n"
        "\tNi = %d, Li = %d, ki = %g, li = %d\n"
        "\tNf = %d, Lf = %d, kf = %g, lf = %d\n",
        lambda, Ni, Li, ki, li, Nf, Lf, kf, lf
    );
    
    if (lambda == 0)
    {
        //
        // r1 > r2
        //
        
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf](double r2) -> double
        {
            // inner integrand
            auto iintegrand = [Ni,Li,kf,lf,r2](double r1) -> double
            {
                return Hydrogen::P(Ni,Li,r1) * special::ric_j(lf,kf*r1) * (1./r1 - 1./r2);
            };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            
            // integrate and check success
            if (not Qi.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Iexc (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }
            
            return Qi.result() * special::ric_j(li,ki*r2) * Hydrogen::P(Nf,Lf,r2);
        };
        
        // outer integrator
        GaussKronrod<decltype(integrand)> Q(integrand);
        
        // integrate and check success
        if (not Q.integrate(0.,Inf))
        {
            throw exception
            (
                "compute_Iexc (outer) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                Ni, Li, ki, li, Nf, Lf, kf, lf, Q.status().c_str(), Q.result()
            );
        }
        
        std::cout << "\tIexc = " << Q.result() << std::endl;
        return Q.result();
    }
    else
    {
        // r1 < r2
        auto integrand1 = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r2) -> double
        {
            // inner integrand
            auto iintegrand1 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double
            {
                return Hydrogen::P(Ni,Li,r1) * special::ric_j(lf,kf*r1) * std::pow(r1/r2,lambda);
            };
            
            // inner integrator
            BesselNodeIntegrator1D<decltype(iintegrand1),GaussKronrod<decltype(iintegrand1)>> Qi1 (iintegrand1, kf, lf);
            
            // integrate and check success
            if (not Qi1.integrate(0.,r2))
            {
                throw exception
                (
                    "compute_Iexc (inner1) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi1.status().c_str(), Qi1.result()
                );
            }
            
            return Qi1.result() / r2 * special::ric_j(li,ki*r2) * Hydrogen::P(Nf,Lf,r2);
        };
        
        BesselNodeIntegrator1D<decltype(integrand1),GaussKronrod<decltype(integrand1)>> Q1 (integrand1, ki, li);
        if (not Q1.integrate(0.,Inf))
        {
            throw exception
            (
                "compute_Iexc (outer1) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q1.status().c_str(), Q1.result()
            );
        }
        
        // r2 < r1
        auto integrand2 = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r1) -> double
        {
            // inner integrand
            auto iintegrand2 = [Nf,Lf,ki,li,r1,lambda](double r2) -> double
            {
                return Hydrogen::P(Nf,Lf,r2) * special::ric_j(li,ki*r2) * std::pow(r2/r1,lambda);
            };
            
            // inner integrator
            BesselNodeIntegrator1D<decltype(iintegrand2),GaussKronrod<decltype(iintegrand2)>> Qi2 (iintegrand2, ki, li);
            
            // integrate and check success
            if (not Qi2.integrate(0.,r1))
            {
                throw exception
                (
                    "compute_Iexc (inner2) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r1, Qi2.status().c_str(), Qi2.result()
                );
            }
            
            return Qi2.result() / r1 * special::ric_j(lf,kf*r1) * Hydrogen::P(Ni,Li,r1);
        };
        
        BesselNodeIntegrator1D<decltype(integrand2),GaussKronrod<decltype(integrand2)>> Q2 (integrand2, kf, lf);
        if (not Q2.integrate(0.,Inf))
        {
            throw exception
            (
                "compute_Iexc (outer2) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q2.status().c_str(), Q2.result()
            );
        }
        
        std::cout << "\tIexc = " << Q1.result() + Q2.result() << std::endl;
        return Q1.result() + Q2.result();
    }
}
