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

#ifndef HEX_DWBA_DWBA2
#define HEX_DWBA_DWBA2

#include "arrays.h"
#include "potential.h"
#include "chebyshev.h"
#include "specf.h"

/**
 * Class members compute contributions to the second order of the distorted
 * wave Born approximation.
 */
class DWBA2
{
public:

    /**
     * Main DWBA2-part function that will compute a scattering amplitude
     * contribution for given angular numbers.
     */
    static void DWBA2_Ln (
        double Ei, int li, int lf, double ki, double kf,
        int Ni, int Nf, int Li, int Lf, int Ln,
        DistortingPotential const & Ui,
        DistortingPotential const & Uf,
        HydrogenFunction    const & psii,
        HydrogenFunction    const & psif,
        DistortedWave       const & chii,
        DistortedWave       const & chif,
        cArray & DD_lf_li_Ln,
        cArray & DE_lf_li_Ln,
        cArray & ED_lf_li_Ln,
        cArray & EE_lf_li_Ln,
        bool & compute_DD,
        bool & compute_DE,
        bool & compute_ED,
        bool & compute_EE
    );
	
	static void DWBA2_energy_driver (
		double Ei, int li, int lf, double ki, double kf, 
		int Ni, int Nf, int Li, int Lf, int Ln,
		int lami, int lamf, 
		HydrogenFunction const & psii, 
		HydrogenFunction const & psif, 
		DistortingPotential const & Ui,
		DistortingPotential const & Uf,
		DistortedWave const & chii,
		DistortedWave const & chif, Complex & cDD
	);
	
	static void DWBA2_En (
		double Ei, int Ni,
		double Kn, int Ln,
		int lami, int lamf,
		HydrogenFunction const & psii,
		HydrogenFunction const & psif,
		HydrogenFunction const & psin,
		DistortingPotential const & Ui,
		DistortingPotential const & Uf,
		DistortingPotential const & Ug,
		DistortedWave const & chii,
		DistortedWave const & chif,
		Complex & DD
	);
	
	/**
	 * Auxiliary class holding information about intermediate integration.
	 * When evaluated, it shall return
	 * \f[
	 *     \int 
	 *         \psi_n(r_1) \left(
	 *             \frac{r_<^\lambda}{r_>^{\lambda+1}} 
	 *           - \frac{\delta_{\lambda 0}}{r_2} 
	 *           - \delta_{\lambda 0} U_\alpha(r_2)
	 *         \right) \psi_\alpha(r_1) \mathrm{d}r_1 \ .
	 * \f]
	 */
	class PhiFunctionDirIntegral : public RadialFunction<double>
	{
	public:
		
		PhiFunctionDirIntegral (
			HydrogenFunction const & psin, 
			int lam, 
			DistortingPotential const & U, 
			HydrogenFunction const & psi
		);
		
		/// evaluated integral
		double operator() (double x2) const;
		
	private:
		
		/// λ
		int Lam;
		
		/// U
		DistortingPotential U;
		
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
	 *         \right) \psi_\alpha(r_1) \mathrm{d}r_1 \ .
	 * \f]
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
		
		/// Chebyshev approximation evaluated at infinity.
		double Cb_inf;
		
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
		
		/// integrand
		PhiFunctionDirIntegral Integrand;
		
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
};

#endif

