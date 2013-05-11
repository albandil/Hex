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
) : Lam(lam), U(U), Diag(psin == psi), 
	Zero(lam == 0 and Diag and (DistortingPotential(psi) == U)),
	Wavenum(psin.getK())
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
// 		if (psin.getK() < 1.)
		{
			// for low energies this is the only way, so we disable the limit
			tryRealChebyshev(psin, psi, -1, Cheb_L_conv, Cheb_mLm1_conv);
		}
// 		else
// 		{
			// for high energies we may back away after finite cblimit
// 			tryRealChebyshev(psin, psi, cblimit, Cheb_L_conv, Cheb_mLm1_conv);
			
			// try to find the expansions using the complex functions
// 			tryComplexChebyshev(psin, psi, cblimit, Cheb_L_conv, Cheb_mLm1_conv);
// 		}
		
		// check success
		if (not Cheb_L_conv or not Cheb_mLm1_conv)
			throw exception ("[PhiFunctionDir] Non-convergent approximation. Try to raise the cblimit.");
		
		// store precomputed expansions to disk
		Cheb_L.coeffs().hdfsave(name1.c_str());
		Cheb_mLm1.coeffs().hdfsave(name2.c_str());
	}
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
	for (int N = 16; (cblimit < 0 or N <= cblimit) and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
	{
		if (not Cheb_L_conv)
		{
			Cheb_L.generate(compact1, N);
			if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
			{
				Cheb_L_conv = true;
				Cheb_L = Cheb_L.integrate(Cheb_L.Integ_Low);
			}
		}
		if (not Cheb_mLm1_conv)
		{
			Cheb_mLm1.generate(compact2, N);
			if ((Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10)) < N)
			{
				Cheb_mLm1_conv = true;
				Cheb_mLm1 = Cheb_mLm1.integrate(Cheb_mLm1.Integ_High);
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
	std::cout << " try complex... " << std::flush;
	
	// expressions to approximate
	auto expr1 = [&](double r2) -> double {
		
		// complex integrand
		auto inte1 = [&](double xi) -> double {
			Complex r1(r2,xi);
			Complex eval = psin.Hplus(r1) * psi(r1) * pow(r1,Lam);
			return eval.real();
		};
		
		// compute the integral
		return -ClenshawCurtis<decltype(inte1),double>(inte1).integrate(0.,Inf);
		
	};
	
	auto expr2 = [&](double r2) -> double {
		
		// complex integrand
		auto inte2 = [&](double xi) -> double {
			Complex r1(r2,xi);
			Complex eval = psin.Hplus(r1) * psi(r1) * pow(r1,-Lam-1);
			return eval.real();
		};
		
		// compute the integral
		return ClenshawCurtis<decltype(inte2),double>(inte2).integrate(0,Inf);
		
	};
	
	// compactifications
	CompactificationR<decltype(expr1),double> compactA(expr1);
	CompactificationR<decltype(expr2),double> compactB(expr2);

	// Coulomb wave prefactor which is not included in psin.Hplus()
	double prefactor = sqrt(M_2_PI)/psin.getK();
	
	// run the evaluation/convergence loop
	for (int N = 16; (cblimit < 0 or N <= cblimit) and not (Cheb_L_conv and Cheb_mLm1_conv); N *= 2)
	{
		std::cout << " " << N << std::flush;
		
		if (not Cheb_L_conv)
		{
			Cheb_L.generate(compactA, N);
			
			std::ostringstream oss;
			oss << "cb_A_" << N;
			Cheb_L.coeffs().hdfsave(oss.str().c_str());
			
			if ((Cheb_L_tail = Cheb_L.tail(1e-10)) < N)
			{
				// multiply all coefficients by missing coefficient sqrt(2/π)/k
				Cheb_L = Chebyshev<double,double>(prefactor * Cheb_L.coeffs());
				Cheb_L_conv = true;
			}
		}
		if (not Cheb_mLm1_conv)
		{
			std::ostringstream oss;
			oss << "cb_B_" << N;
			Cheb_L.coeffs().hdfsave(oss.str().c_str());
			
			Cheb_mLm1.generate(compactB, N);
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
	
	return pow(r,Lam) * Cheb_mLm1.clenshaw(x, Cheb_mLm1_tail)
		+ pow(r, -Lam-1) * Cheb_L.clenshaw(x, Cheb_L_tail);
}
