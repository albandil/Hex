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

#ifndef HEX_DWBA_MULTI
#define HEX_DWBA_MULTI

#include "arrays.h"
#include "hydrogen.h"
#include "potential.h"
#include "chebyshev.h"
#include "specf.h"
	
/**
 * Auxiliary class holding information about intermediate integration.
 * When evaluated, it shall return
 * \f[
 *     \int 
 *         \psi_n(r_1) \left(
 *             \frac{r_<^\lambda}{r_>^{\lambda+1}} 
 *           - \frac{\delta_{\lambda 0}}{r_2} 
 *           - \delta_{\lambda 0} U_\alpha(r_2)
 *         \right) \psi_\alpha(r_1) \mathrm{d}r_1 
 * \f]
 * as a result of Chebyshev approximation of \ref PhiFunctionDirIntegral .
 */
class PhiFunctionDir : public RadialFunction<double>
{
public:
	
	PhiFunctionDir (
		HydrogenFunction const & psin, 
		int lam, 
		DistortingPotential const & U, 
		HydrogenFunction const & psi
	);
	
	/// Evaluate the function.
	double operator() (double x) const;
	
private:
	
	int Lam;
	DistortingPotential U;
	
	/// whether psin == psi
	bool Diag;
	
	/// whether the integral is identical zero for all "x2"
	bool Zero;
	
	/// Chebyshev approximations of the integrand \f$ \psi_n(r) r^\lambda \psi_i(r) \f$.
	Chebyshev<double,double> Cheb_L;
	
	/// Chebyshev approximations of the integrand \f$ \psi_n(r) r^{-\lambda-1} \psi_i(r) \f$.
	Chebyshev<double,double> Cheb_mLm1;
	
	int Cheb_mLm1_tail;
	int Cheb_L_tail;
	
	double Cheb_mLm1_inf;
	double Cheb_L_zero;
};

#endif
