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

class PhiFunction : public RadialFunction<double>
{
	public:
		
		virtual double operator() (double x) const = 0;
		virtual bool isZero() const = 0;
		virtual double k() const = 0;
};

/**
 * \brief Multipole function of the "direct" type.
 * 
 * Auxiliary class holding information about intermediate integration.
 * When evaluated, it shall return
 * \f[
 *     \int 
 *         \psi_n(r_1) \left(
 *             \frac{r_<^\lambda}{r_>^{\lambda+1}} 
 *           - \frac{\delta_{\lambda 0}}{r_2} 
 *           - \delta_{\lambda 0} U_\alpha(r_2)
 *         \right) \psi_\alpha(r_1) \mathrm{d}r_1  \ .
 * \f]
 * 
 * The integral is split into two integrals with a fixed integrand (no \f$ r_< \f$,
 * \f$ r_> \f$ anymore). The constructor tries to approximate the integrands
 * using a standard Chebyshev approximation. If successfull, the expansions are
 * trivially integrated and stored for use in evaluation by the operator().
 * 
 * If the maximal given number of allowed Chebyshev nodes to use is hit, then
 * a different method is used. Again, a Chebyshev approximation is sought, but
 * now of the whole integral as a function of \f$ r_2 \f$. The integral is done
 * by integration in complex plane, along a contour, where the oscillating free
 * intermediate state quickly damps.
 * 
 * \param psin Intermediate hydrogen function (bound or free).
 * \param lam Multipole moment \f$ \lamda \ge 0 \f$.
 * \param U Distorting potential.
 * \param psi Initial (or final) atomic bound state.
 * \param cblimit Maximal number of Chebyshev evaluation points used. If the
 *     routine doesn't find a satisfactory approximation even after evaluating
 *     the function at 'cblimit' Chebyshev roots, the integrand is assumed
 *     to be highly oscillatory and an alternative method is used. Minus one disables
 *     the limit.
 */
class PhiFunctionDir : public PhiFunction
{
public:
	
	PhiFunctionDir (
		HydrogenFunction const & psin, 
		int lam, 
		DistortingPotential const & U, 
		HydrogenFunction const & psi,
		int cblimit = 1024
	);
	
	/// Evaluate the function.
	double operator() (double x) const;
	
	/// Whether the function is identical zero.
	bool isZero() const { return Zero; }
	
	/// Return effective wavenumber.
	double k() const { return Wavenum; }
	
private:
	
	/// Compute Chebyshev expansion coefficients by evaluating the integrands.
	void tryFullRealChebyshev (
		int cblimit,
		bool & Cheb_L_conv,
		bool & Cheb_mLm1_conv
	);
	
	void tryFrontRealChebyshev (
		int cblimit,
		bool & Cheb_L_conv,
		bool & Cheb_mLm1_conv
	);
	
	/// Compute Chebyshev expansion coefficients by computing the complex integral.
	void tryFullComplexChebyshev (
		int cblimit,
		bool & Cheb_L_conv,
		bool & Cheb_mLm1_conv
	);
	
	HydrogenFunction psin;
	HydrogenFunction psi;
	
	int Lam;
	DistortingPotential U;
	
	/// whether psin == psi
	bool Diag;
	
	/// whether the integral is identical zero for all "x2"
	bool Zero;
	
	/// Wavenumber of the intermediate state.
	double Wavenum;
	
	/// Use frontal approximation.
	bool UseFront;
	
	/// Chebyshev approximations of the integrand \f$ \psi_n(r) r^\lambda \psi_i(r) \f$.
	Chebyshev<double,double> Cheb_L;
	
	/// Chebyshev approximations of the integrand \f$ \psi_n(r) r^{-\lambda-1} \psi_i(r) \f$.
	Chebyshev<double,double> Cheb_mLm1;
	
	Chebyshev<double,double> ACoeff, BCoeff;
	int ATail, BTail;
	
	int Cheb_mLm1_tail;
	int Cheb_L_tail;
};

#endif
