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
 * as a result from direct integration.
 */
class PhiFunctionDirIntegral : public RadialFunction<double>
{
public:
	
	PhiFunctionDirIntegral (
		HydrogenFunction const & psin, 
		int lam,
		HydrogenFunction const & psi
	);
	
	/// evaluated integral
	double operator() (double x2) const;
	
private:
	
	/// λ
	int Lam;
	
	/// psi
	HydrogenFunction Psi;
	
	/// psin
	HydrogenFunction Psin;
};
	
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
	
	/// Compose filename of the HDF storage file.
	static std::string name(HydrogenFunction const & psin, int lam, HydrogenFunction const & psi);
	
	/// Load from a HDF file.
	bool load(std::string filename);
	
	/// Save to a HDF file.
	void save(std::string filename) const;
	
	/// Truncation index of the Chebyshev approximation.
	int Tail;
	
	/// λ (a constructor parameter)
	int Lam;
	
	/// U (a constructor parameter)
	DistortingPotential U;
	
	/// whether psin == psi
	bool Diag;
	
	/// whether the integral is identical zero for all "x2"
	bool Zero;
	
	/// integranl
	PhiFunctionDirIntegral Integral;
	
	/**
		* Compactification of the function
		* \f[
		*     \Phi(r_2) = \int_0^\infty \psi_n(r_1) \left(
		*          \frac{r_<^\lambda}{r_>^{\lambda+1}}
		*        - \delta_{\lambda 0} \frac{1}{r_2}
		*        - \delta_{\lambda 0} U_\alpha(r_2)
		*     \right) \psi_\alpha(r_1) \mathrm{d}r_1 \ .
		* \f]
		*/
	CompactificationR<PhiFunctionDirIntegral,double> CompactIntegral;
	
	/// Chebyshev approximation of the compactified function \f$ \Phi(r_2) \f$.
	Chebyshev<double,double> CompactIntegralCb;
};

#endif
