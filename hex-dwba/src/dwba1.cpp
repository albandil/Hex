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

#include "arrays.h"
#include "complex.h"
#include "potential.h"
#include "dwba1.h"
#include "hydrogen.h"
#include "gausskronrod.h"
#include "specf.h"

namespace DWBA1
{

Complex computeDirect1e(DistortingPotential const& U, int l, double k)
{
	// get distorted wave
	DistortedWave chi_kl = U.getDistortedWave(k,l);
	
	// set up the integrands
	auto integrand1 = [ & ](double r) -> double {
		return chi_kl(r) * U(r) * ric_j(l,k*r);
	};
	auto integrand2 = [ & ](double r) -> double {
		return chi_kl(r) * U(r) * chi_kl(r);
	};
	
	// integrate
	GaussKronrod<decltype(integrand1)> Q1(integrand1);
	GaussKronrod<decltype(integrand2)> Q2(integrand2);
	Q1.integrate(0., Inf);
	Q2.integrate(0., Inf);
	
	// return the result
	return pow(4*M_PI,2) * sqrt((2*l+1)/(4*M_PI)) / (k*k) * 
		(Q1.result() - Q2.result() * chi_kl.getPhasef()) * chi_kl.getPhasef();
}

Complex computeExchange1e (
	DistortingPotential const& U,
	int Ni, int Li, double ki,
	int Nf, int Lf, double kf
) {
	
	// get distorted waves
	DistortedWave chif = U.getDistortedWave(kf,Li);
	DistortedWave chii = U.getDistortedWave(ki,Lf);
	
	// set up the integrands
	auto integrand1 = [ & ](double r) -> double {
		return chif(r) * U(r) * Hydrogen::evalBoundState(Ni,Li,r);
	};
	auto integrand2 = [ & ](double r) -> double {
		return Hydrogen::evalBoundState(Nf,Lf,r) * chii(r);
	};
	
	// integrate
	GaussKronrod<decltype(integrand1)> Q1(integrand1);
	GaussKronrod<decltype(integrand2)> Q2(integrand2);
	Q1.integrate(0., Inf);
	Q2.integrate(0., Inf);
	
	// compute phase factor
	double phase = chif.getPhase() + chii.getPhase();
	Complex phasef(cos(phase),sin(phase));
	
	// return the result
	return -pow(4*M_PI,2) * phasef * pow(Complex(0.,1.), Lf-Li) * 
		sqrt((2*Lf+1)/(4*M_PI)) * Q1.result() * Q2.result() / (ki*kf);
}

Complex computeDirect2e (
	const DistortingPotential& U, int lambda,
	int Nf, int Lf, double kf, int lf,
	int Ni, int Li, double ki, int li
) {
	// get distorted waves
	DistortedWave chif = U.getDistortedWave(kf,lf);
	DistortedWave chii = U.getDistortedWave(ki,li);
	
	// set up the outer integrand
	auto outer_integrand = [ & ](double r2) -> double {
		
		// set up the inner integrands
		auto inner_integrand_1 = [ & ](double r1) -> double {
			double psi_f = Hydrogen::evalBoundState(Nf,Lf,r1);
			double psi_i = Hydrogen::evalBoundState(Ni,Li,r1);
			double multipole;
			if (lambda == 0)
				multipole = 1/r1 - 1/r2;
			else
				multipole = pow(r2/r1,lambda)/r1;
			return psi_f * multipole * psi_i;
		};
		auto inner_integrand_2 = [ & ](double r1) -> double {
			double psi_f = Hydrogen::evalBoundState(Nf,Lf,r1);
			double psi_i = Hydrogen::evalBoundState(Ni,Li,r1);
			double multipole;
			if (lambda == 0)
				return 0;
			else
				multipole = pow(r1/r2,lambda)/r2;
			return psi_f * multipole * psi_i;
		};
		
		// integrate
		GaussKronrod<decltype(inner_integrand_1)> Q1(inner_integrand_1);
		GaussKronrod<decltype(inner_integrand_2)> Q2(inner_integrand_2);
		Q1.integrate(r2, Inf);
		Q2.integrate(0, r2);
		
		// evaluate distorted waves
		double chi_i = chii(r2);
		double chi_f = chif(r2);
		
		// return the result
		return chi_i * chi_f * (Q1.result() + Q2.result());
	};
	
	// integrate
	GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
	Q.integrate(0., Inf);
	
	// compute phase factor
	double phase = chii.getPhase() + chif.getPhase();
	Complex phasef(cos(phase),sin(phase));
	
	// return the result
	return pow(4*M_PI,3) * phasef * pow(Complex(0.,1.),li-lf) * 
		sqrt((2*li+1)/(4*M_PI)) / (2.*lambda+1.) * Q.result() / (ki*kf);
}

Complex computeExchange2e (
	const DistortingPotential& U, int lambda,
	int Nf, int Lf, double kf, int lf,
	int Ni, int Li, double ki, int li
) {
	// get distorted waves
	DistortedWave chif = U.getDistortedWave(kf,lf);
	DistortedWave chii = U.getDistortedWave(ki,li);
	
	// set up the outer integrand
	auto outer_integrand = [ & ](double r2) -> double {
		
		// set up the inner integrands
		auto inner_integrand_1 = [ & ](double r1) -> double {
			double psi_f = Hydrogen::evalBoundState(Nf,Lf,r1);
			double chi_i = chii(r1);
			double multipole;
			if (lambda == 0)
				return 0.;
			else
				multipole = pow(r2/r1,lambda)/r1;
			return psi_f * multipole * chi_i;
		};
		auto inner_integrand_2 = [ & ](double r1) -> double {
			double psi_f = Hydrogen::evalBoundState(Nf,Lf,r1);
			double chi_i = chii(r1);
			double multipole;
			if (lambda == 0)
				multipole = 1/r2 - 1/r1;
			else
				multipole = pow(r1/r2,lambda)/r2;
			return psi_f * multipole * chi_i;
		};
		
		// integrate
		GaussKronrod<decltype(inner_integrand_1)> Q1(inner_integrand_1);
		GaussKronrod<decltype(inner_integrand_2)> Q2(inner_integrand_2);
		Q1.integrate(r2, Inf);
		Q2.integrate(0, r2);
		
		// evaluate distorted waves
		double psi_i = Hydrogen::evalBoundState(Ni,Li,r2);
		double chi_f = chif(r2);
		
		// return the result
		return psi_i * chi_f * (Q1.result() + Q2.result());
	};
	
	// integrate
	GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
	Q.integrate(0., Inf);
	
	// compute phase factor
	double phase = chii.getPhase() + chif.getPhase();
	Complex phasef(cos(phase),sin(phase));
	
	// return the result
	return pow(4*M_PI,3) * phasef * pow(Complex(0.,1.),li-lf) * 
		sqrt((2*li+1)/(4*M_PI)) / (2.*lambda+1.) * Q.result() / (ki*kf);
}

}; // end of namespace DWBA1
