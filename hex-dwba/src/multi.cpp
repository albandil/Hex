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
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
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
) : Lam(lam), U(U), Diag(psin == psi), 
	Zero(lam == 0 and Diag and (DistortingPotential(psi) == U))
{
	if (Zero)
		return;
	
	std::ostringstream oss;
	oss << "phidir-" << lam << "-"
	    << psi.getK()  << "-" << psi.getN()  << "-" << psi.getL()  << "-"
	    << psin.getK() << "-" << psin.getN() << "-" << psin.getL() << "~";
	std::string name1 = oss.str() + "1.chb";
	std::string name2 = oss.str() + "2.chb";
	
	rArray a, b;
	if (a.hdfload(name1.c_str()) and b.hdfload(name2.c_str()))
	{
		Cheb_L = Chebyshev<double,double>(a);
		Cheb_mLm1 = Chebyshev<double,double>(b);
		
		Cheb_L_tail = Cheb_L.tail(1e-10);
		Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10);
	}
	else
	{
		// auxiliary variables indicating success when approximating the functions
		bool Cheb_L_conv = false, Cheb_mLm1_conv = false;
		
		// try to find the expansions using the real functions
		tryRealChebyshev(psin, psi, cblimit, Cheb_L_conv, Cheb_mLm1_conv);
		
		// try to find the expansions using the complex functions
		tryComplexChebyshev(psin, psi, cblimit, Cheb_L_conv, Cheb_mLm1_conv);
		
		// check success
		if (not Cheb_L_conv or not Cheb_mLm1_conv)
			throw exception ("[PhiFunctionDir] Non-convergent approximation. Try to raise the cblimit.");
		
		// store precomputed expansions to disk
		Cheb_L.coeffs().hdfsave(name1.c_str());
		Cheb_mLm1.coeffs().hdfsave(name2.c_str());
	}
	
	Cheb_mLm1_inf = Cheb_mLm1.clenshaw(1, Cheb_mLm1_tail);
	Cheb_L_zero = Cheb_L.clenshaw(-1, Cheb_L_tail);
}

void PhiFunctionDir::tryRealChebyshev (
	HydrogenFunction const & psin, 
	HydrogenFunction const & psi,
	int cblimit,
	bool & Cheb_L_conv,
	bool & Cheb_mLm1_conv
){
	// integrands
	auto inte1 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,Lam); };
	auto inte2 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,-Lam-1); };
	
	// integrand compactifications
	CompactIntegrand<decltype(inte1),double> compact1(inte1);
	CompactIntegrand<decltype(inte2),double> compact2(inte2);
	
	// run the evaluation/convergence loop
	for (int N = 16; N <= cblimit and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
	{
		if (not Cheb_L_conv)
		{
			Cheb_L.generate(compact1, N);
			if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
			{
				Cheb_L_conv = true;
				Cheb_L = Cheb_L.integrate();
			}
		}
		if (not Cheb_mLm1_conv)
		{
			Cheb_mLm1.generate(compact2, N);
			if ((Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10)) < N)
			{
				Cheb_mLm1_conv = true;
				Cheb_mLm1 = Cheb_mLm1.integrate();
			}
		}
	}
}

void PhiFunctionDir::tryComplexChebyshev (
	HydrogenFunction const & psin, 
	HydrogenFunction const & psi,
	int cblimit,
	bool & Cheb_L_conv,
	bool & Cheb_mLm1_conv
){
	// expressions to approximate
	auto expr1 = [&](double r2) -> double {
		
		// complex integrand
		auto inte1 = [&](double xi) -> double {
			Complex r1(r2,xi);
			return (psin.Hplus(r1)*psi(r1)*pow(r1,Lam)).real();
		};
		
		// integrand compactification
		CompactIntegrand<decltype(inte1),double> compact1(inte1);
		
		// compute the integral
		ClenshawCurtis<decltype(compact1),double> ccint(compact1);
		return ccint.integrate (
			compact1.scale(0.),
			compact1.scale(r2)
		);
		
	};
		
	auto expr2 = [&](double r2) -> double {
		
		// complex integrand
		auto inte2 = [&](double xi) -> double {
			Complex r1(r2,xi);
			return (psin.Hplus(r1)*psi(r1)*pow(r1,-Lam-1)).real();
		};
		
		// integrand compactification
		CompactIntegrand<decltype(inte2),double> compact2(inte2);
		
		// compute the integral
		ClenshawCurtis<decltype(compact2),double> ccint(compact2);
		return ccint.integrate (
			compact2.scale(r2),
			compact2.scale(Inf)
		);
		
	};
	
	// Coulomb wave prefactor which is not included in psin.Hplus()
	double prefactor = sqrt(M_2_PI)/psin.getK();
	
	// run the evaluation/convergence loop
	for (int N = 16; N <= cblimit and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
	{
		if (not Cheb_L_conv)
		{
			Cheb_L.generate(expr1, N);
			if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
			{
				// multiply all coefficients by missing coefficient sqrt(2/π)/k
				Cheb_L = Chebyshev<double,double>(prefactor * Cheb_L.coeffs());
				Cheb_L_conv = true;
			}
		}
		if (not Cheb_mLm1_conv)
		{
			Cheb_mLm1.generate(expr2, N);
			if ((Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10)) < N)
			{
				// multiply all coefficients by missing coefficient sqrt(2/π)/k
				Cheb_mLm1 = Chebyshev<double,double>(prefactor * Cheb_mLm1.coeffs());
				Cheb_mLm1_conv = true;
			}
		}
	}
}

double PhiFunctionDir::operator() (double r) const
{
	double x = (r - 1) / (r + 1);
	
	if (Zero)
		return 0.;
	
	if (r == 0)
		return 0.;	// FIXME
	
	return pow(r,Lam) * (Cheb_mLm1_inf - Cheb_mLm1.clenshaw(x, Cheb_mLm1_tail))
		+ pow(r, -Lam-1) * (Cheb_L.clenshaw(x, Cheb_L_tail) - Cheb_L_zero);
}