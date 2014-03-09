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

#include <deque>
#include <iostream>
#include <vector>

#include "arrays.h"
#include "gausskronrod.h"
#include "radial.h"
#include "romberg.h"
#include "specf.h"
#include "symbolic.h"

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
                
                std::cout << "integrating node : " << rmin << " to " << rmax << " : " << Q_.result() << std::endl;
                
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

template <class Functor, class Integrator>
class BesselNodeIntegrator2D
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
        
        typedef struct { int x,y; } tNode;
        
        BesselNodeIntegrator2D (Functor f, double k, int l)
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
        
        /**
         * @brief Integrate on @f$ [x_1,x_2] \times [y_1,y_2] @f$.
         */
        bool integrate (double x1, double y1, double x2, double y2)
        {
            std::deque<tNode> nodes;
            std::vector<double> xroots = { 0. };
            std::vector<double> yroots = { 0. };
            
            // add initial integration parcel
            nodes.push_back({0,0});
            
            // value of the integral
            double integral = 0.;
            
            // while there are any parcels for integration
            for (int i = 0; not nodes.empty(); i++)
            {
                // get front parcel
                tNode node = nodes.pop_front();
                
                // update root cache (if necessary)
                if (node.x + 1 == xroots.size())
                {
                    // compute new Bessel function node
                    // TODO
                }
                if (node.y + 1 == yroots.size())
                {
                    // compute new Bessel function node
                    // TODO
                }
                
                // get parcel integration nodes
                double xmin = std::max (xroots[node.x],     x1);
                double xmax = std::min (xroots[node.x + 1], x2);
                double ymin = std::max (yroots[node.y],     y1);
                double ymax = std::min (yroots[node.y + 1], y2);
                
                // check size of the parcel (skip if zero-sized)
                if (xmin >= xmax or ymin >= ymax)
                    continue;
                
                // integrate on the parcel [xmin,xmax] × [ymin,ymax]
                double contrib = Q_.integrate (xmin, xmax, ymin, ymax);
                integral += contrib;
                
                // decide whether the contribution from this parcel is significant
                if (std::abs(contrib) < epsabs_ or std::abs(contrib) < epsrel_ * std::abs(integral))
                {
                    // add surrounding nodes to the queue
                    nodes.push_back({node.x+1, node.y  });
                    nodes.push_back({node.x,   node.y+1});
                    nodes.push_back({node.x+1, node.y+1});
                }
                
                // check iteration limit
                if (i == limit_)
                {
                    ok_ = false;
                    status_ = "iteration limit reached";
                    result_ = integral;
                    return ok_;
                }
            }
            
            ok_ = true;
            status_ = "";
            result_ = integral;
            return ok_;
        }
};

double inner_Idir (int lambda, int Ni, int Li, int Nf, int Lf, double r, double sumtol = 1e-10)
{
    // compute normaliazation factors of the generalized Laguerre polynomials
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
    
    // inner integral
    double integral = 0.;
    
    // factor in the argument of the exponential
    double c = 1./Ni + 1./Nf;
    
    // integrate term by term
    for (symbolic::term p : PP)
    {
        // check if the integral is divergent
        if (p.a < lambda + 1)
        {
            throw exception
            (
                "[inner_Idir] Integral is divergent, a < λ + 1, for a = %d, λ = %d.",
                p.a, lambda
            );
        }
        
        
        // compute high integral
        double cr_to_k = 1;
        double k_fac = 1;
        double suma1 = 1;
        for (int k = 1; k <= p.a - lambda - 1; k++)
        {
            cr_to_k *= c * r;
            k_fac *= k;
            
            suma1 += cr_to_k / k_fac;
        }
        suma1 *= gsl_sf_fact(p.a-lambda-1) / gsl_sf_pow_int(c,p.a-lambda);
        
        // compute low integral
        cr_to_k = gsl_pow_int(r,p.a-lambda);
        double poch = p.a + lambda + 1;
        double suma2 = 1/poch;
        for /*ever*/ (int k = 1; ; k++)
        {
            cr_to_k *= c * r;
            poch *= p.a + lambda + 1 + k;
            
            double contrib = cr_to_k / poch;
            suma2 += contrib;
            
            if (std::abs(contrib) < sumtol * std::abs(suma1 + suma2))
                break;
        }
        
        // update integral
        integral += (suma1 + suma2) * symbolic::double_approx(p.kr);
    }
    
    return Normi * Normf * integral * std::pow(r, lambda) * std::exp(-c*r);
}

double compute_Idir (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
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
            auto iintegrand = [Ni,Li,Nf,Lf,r2](double r1) -> double { return hydro_P(Ni,Li,r1) * hydro_P(Nf,Lf,r1) * (1./r1 - 1./r2); };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand)> Qi(iintegrand);
            
            // integrate and check success
            if (not Qi.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Idir (inner) failed for λ=0, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Qi.status().c_str(), Qi.result()
                );
            }
            
            return Qi.result() * ric_j(li,ki*r2) * ric_j(lf,kf*r2);
        };
        
        // which Bessel function oscillates slowlier ?
        int    l = (ki < kf ? li : lf);
        double k = (ki < kf ? ki : kf);
        
        // outer integrator
        BesselNodeIntegrator1D<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.integrate(0, Inf);
        
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
        
        std::cout << PP << std::endl;
        std::cout << Normi << std::endl;
        std::cout << Normf << std::endl;
        
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
                    int_high = gsl_sf_pow_int(r,lambda) * gsl_sf_pow_int(c,lambda-p.a) * res.val;
                
                // compute the low integral
                int err_low = gsl_sf_gamma_inc_P_e (p.a + lambda - 1, c * r, &res);
                double scale = gsl_sf_gamma (p.a + lambda - 1);
                if (err_low != GSL_SUCCESS)
                {
                    throw exception
                    (
                        "Unable to evaluate scaled complementary incomplete gamma-function P(%d,%g) - %s.",
                        p.a + lambda - 1, c * r, gsl_strerror(err_low)
                    );
                }
                double int_low = gsl_sf_pow_int(r,-lambda-1) * gsl_sf_pow_int(c,-p.a-lambda-1) * res.val * scale;
                
                // sum both contributions
                integral = (int_low + int_high) * symbolic::double_approx(p.kr);
                
//                 std::cout << "p: " << p << " : " << (int_low + int_high) * symbolic::double_approx(p.kr) << std::endl;
            }
            
            return Normi * Normf * integral * ric_j(li,ki*r) * ric_j(lf,kf*r);
        };
        
        // which Bessel function oscillates slowlier ?
        int    l = (ki < kf ? li : lf);
        double k = (ki < kf ? ki : kf);
        
        // outer integrator
        BesselNodeIntegrator1D<decltype(integrand),GaussKronrod<decltype(integrand)>> R(integrand, k, l);
        R.integrate (0, Inf);
        std::cout << "Integral: " << R.result() << std::endl;
        
        return R.result();
    }
}

double compute_Iexc (int li, int lf, int lambda, int Ni, int Li, double ki, int Nf, int Lf, double kf)
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
            auto iintegrand = [Ni,Li,kf,lf,r2](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * (1./r1 - 1./r2); };
            
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
            
            return Qi.result() * ric_j(li,ki*r2) * hydro_P(Nf,Lf,r2);
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
        
        return Q.result();
    }
    else
    {
        /*
        auto integrand = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r2) -> double
        {
            // inner integrand
            auto iintegrand1 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * std::pow(r1/r2,lambda); };
            auto iintegrand2 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * std::pow(r2/r1,lambda+1); };
            
            // inner integrator
            GaussKronrod<decltype(iintegrand1)> Q1 (iintegrand1);
            GaussKronrod<decltype(iintegrand2)> Q2 (iintegrand2);
            
            // integrate and check success
            if (not Q1.integrate(0.,r2))
            {
                throw exception
                (
                    "compute_Iexc (inner1) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Q1.status().c_str(), Q1.result()
                );
            }
            if (not Q2.integrate(r2,Inf))
            {
                throw exception
                (
                    "compute_Iexc (inner2) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d, r2=%g (\"%s\").\n\tresult = %g\n",
                    lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, r2, Q2.status().c_str(), Q2.result()
                );
            }
            
            return (Q1.result() + Q2.result()) / r2 * ric_j(li,ki*r2) * hydro_P(Nf,Lf,r2);
        };
        
        // outer integrator
        GaussKronrod<decltype(integrand)> Q(integrand);
        
        // integrate and check success
        if (not Q.integrate(0.,Inf))
        {
            throw exception
            (
                "compute_Iexc (outer) failed for λ=%d, Ni=%d, Li=%d, ki=%g, li=%d, Nf=%d, Lf=%d, kf=%g, lf=%d (\"%s\").\n\tresult = %g",
                lambda, Ni, Li, ki, li, Nf, Lf, kf, lf, Q.status().c_str(), Q.result()
            );
        }
        
        return Q.result();
        */
        
        // r1 < r2
        auto integrand1 = [Ni,Li,Nf,Lf,li,ki,lf,kf,lambda](double r2) -> double
        {
            // inner integrand
            auto iintegrand1 = [Ni,Li,kf,lf,r2,lambda](double r1) -> double { return hydro_P(Ni,Li,r1) * ric_j(lf,kf*r1) * std::pow(r1/r2,lambda); };
            
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
            
            return Qi1.result() / r2 * ric_j(li,ki*r2) * hydro_P(Nf,Lf,r2);
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
            auto iintegrand2 = [Nf,Lf,ki,li,r1,lambda](double r2) -> double { return hydro_P(Nf,Lf,r2) * ric_j(li,ki*r2) * std::pow(r2/r1,lambda); };
            
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
            
            return Qi2.result() / r1 * ric_j(lf,kf*r1) * hydro_P(Ni,Li,r1);
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
        
        return Q1.result() + Q2.result();
    }
}
