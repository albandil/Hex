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
#include "hydrogen.h"
#include "multi.h"

PhiFunctionDir::PhiFunctionDir (
	HydrogenFunction const & psin, 
	int lam, 
	DistortingPotential const & U, 
	HydrogenFunction const & psi
) : Lam(lam), U(U), Diag(psin == psi), 
	Zero(lam == 0 and Diag and (DistortingPotential(psi) == U))
{
	if (Zero)
		return;
	
	std::ostringstream oss;
	oss << "phidir-" << lam << "-"
	    << psi.getK()  << "-" << psi.getN()  << "-" << psi.getL()  << "-"
	    << psin.getK() << "-" << psin.getN() << "-" << psin.getL() << "~";
	std::string name1 = oss.str() + "1.hdf";
	std::string name2 = oss.str() + "2.hdf";
	
	rArray a, b;
	if (load_array(a, name1.c_str()) and load_array(b, name2.c_str()))
	{
		Cheb_L = Chebyshev<double,double>(a);
		Cheb_mLm1 = Chebyshev<double,double>(b);
		
		Cheb_L_tail = Cheb_L.tail(1e-10);
		Cheb_mLm1_tail = Cheb_mLm1.tail(1e-10);
	}
	else
	{
		bool Cheb_L_conv = false;
		bool Cheb_mLm1_conv = false;
		
		auto inte1 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,Lam); };
		auto inte2 = [&](double x) -> double { return psin(x)*psi(x)*pow(x,-Lam-1); };
		
		CompactIntegrand<decltype(inte1),double> compact1(inte1);
		CompactIntegrand<decltype(inte2),double> compact2(inte2);
		
		for (int N = 16; not Cheb_L_conv or not Cheb_mLm1_conv; N *= 2)
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
		
		save_array(Cheb_L.coeffs(), name1.c_str());
		save_array(Cheb_mLm1.coeffs(), name2.c_str());
	}
	
	Cheb_mLm1_inf = Cheb_mLm1.clenshaw(1, Cheb_mLm1_tail);
	Cheb_L_zero = Cheb_L.clenshaw(-1, Cheb_L_tail);
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
