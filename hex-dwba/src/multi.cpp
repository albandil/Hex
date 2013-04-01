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

#include "multi.h"

PhiFunctionDirIntegral::PhiFunctionDirIntegral (
	HydrogenFunction const & psin, 
	int lam,
	HydrogenFunction const & psi
) : Lam(lam), Psi(psi), Psin(psin) {}

double PhiFunctionDirIntegral::operator()(double x2) const
{
// 	std::cout << "[PhiFunctionDirIntegral::operator()] x2 = " << x2 << "\n";
	if (not finite(x2))
		return 0.;
	
	//
	// finite integrand
	//
		
	auto integrand1 = [&](double x1) -> double {
		if (x1 == 0. or x2 == 0. or not finite(x1))
			return 0.;
		return Psin(x1) * pow(x1/x2, Lam) * Psi(x1);
	};
	
	// compactification
	CompactIntegrand<decltype(integrand1),double> R1(integrand1, 0., Inf, false, 1.0);
	
	// integration system
	ClenshawCurtis<decltype(R1),double> Q1(R1);
	Q1.setEps(1e-10);	// be precise
// 	Q1.setVerbose(true);
	
	// integrate
// 	std::cout << "i1 = ?\n";
	double i1 = (x2 == 0) ? 0. : Q1.integrate(R1.scale(0.), R1.scale(x2)) / x2;
// 	std::cout << "i1 = " << i1 << "\n";
	
	//
	// infinite integrand:
	//
	
	auto integrand2 = [&](double x1) -> double {
		if (x1 == 0. or not finite(x1))
			return 0.;
		return Psin(x1) * pow(x2/x1, Lam) * Psi(x1) / x1;
	};
	
	// compactification
	CompactIntegrand<decltype(integrand2),double> R2(integrand2, 0., Inf, false, 1.0);
	
	// integration system
	ClenshawCurtis<decltype(R2),double> Q2(R2);
	Q2.setEps(1e-10);
	
	// integrate
// 	std::cout << "i2 = ?\n";
	double i2 = Q2.integrate(R2.scale(x2), R2.scale(Inf));
// 	std::cout << "i2 = " << i2 << "\n";
	
	// sum the two integrals
	return i1 + i2;
};

PhiFunctionDir::PhiFunctionDir (
	HydrogenFunction const & psin, 
	int lam, 
	DistortingPotential const & U, 
	HydrogenFunction const & psi
) : Lam(lam), U(U), Diag(psin == psi), 
	Zero(lam == 0 and Diag and (DistortingPotential(psi) == U)),
	Integral(psin,lam,psi),
	CompactIntegral(Integral,0.,false)
{
	if (Zero)
		return;
	
// 	for (int ix = 0; ix <= 1000; ix++)
// 	{
// 		double x = ix * 0.01;
// 		std::cout << x << "\t" 
// 		          << Integral(x) << "\t"
// 		          << CompactIntegral(CompactIntegral.scale(x)) << std::endl;
// 	}
	
	// try to load PhiFunctionDir from a HDF file
	if (not load(name(psin, lam, psi)))
	{	
		// convergence loop
		for (int N = 16; ; N *= 2)
		{
			// construct a Chebyshev approximation of the compactified function
			CompactIntegralCb = Chebyshev<double,double> (CompactIntegral, N);
			
			// check convergence
			if (CompactIntegralCb.tail(1e-10) != N)
				break;
			
			// non-convergent cases need to be done in some other way
			if (N == 1024)
				/* TODO */;
		}
		
		// save PhiFunctionDir to a HDF file
		save(name(psin, lam, psi));
	}
	
	// get optimal truncation index
	Tail = CompactIntegralCb.tail(1e-10);
}

double PhiFunctionDir::operator() (double x) const
{
	if (Zero)
		return 0.;
	
	if (Lam == 0 and Diag)
		return CompactIntegralCb.clenshaw(CompactIntegral.scale(x), Tail) - U.plusMonopole(x);
	
	return CompactIntegralCb.clenshaw(CompactIntegral.scale(x), Tail);
}

std::string PhiFunctionDir::name(HydrogenFunction const & psin, int lam, HydrogenFunction const & psi)
{
	// compose the filename
	//   "phidir-<psin>-<lam>-<psi>.hdf"
	std::ostringstream oss;
	oss << "phidir-"
	    << psin.getN() << "-" << psin.getK() << "-" << psin.getL() << "-"
		<< lam << "-"
		<< psi.getN()  << "-" << psi.getK()  << "-" << psi.getL()
		<< ".arr";
	return oss.str();
}

bool PhiFunctionDir::load(std::string filename)
{
	rArray coeffs;
	
	// load coefficients
	if (load_array(coeffs, filename.c_str()))
	{
		CompactIntegralCb = Chebyshev<double,double> (coeffs);
		return true;
	}
	
	return false;
}

void PhiFunctionDir::save(std::string filename) const
{
	save_array(CompactIntegralCb.coeffs(), filename.c_str());
}
