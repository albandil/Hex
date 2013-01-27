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

#include <cassert>
#include <map>
#include <tuple>

#include "complex.h"
#include "potential.h"
#include "chebyshev.h"
#include "specf.h"
#include "integrate.h"

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
		int Ni, int Nf, int Li, int Lf, int L1, int L2,
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
	
	/**
	 * Class holding an intermediate result from inner integration.
	 */
	class PhiFunction : public RadialFunction<double>
	{
	public:
		
		PhiFunction() {}
		
		PhiFunction (
			double (PhiFunction::* integrator) (double) const,
			int lambda,
			HydrogenFunction const & psi_a,
			SturmianFunction const & S_u
		) {
			// make sure we won't get stuck in cyclical function call
			assert(integrator != &PhiFunction::operator());
			
			// store pointer-to-member to use wen evaluating
			Eval = integrator;
			
			// the integrand we want to approximate by Chebyshev polynomials
			auto integrand = [&](double x) -> double {
				return pow(x,lambda) * psi_a(x) * S_u(x);
			};
			
			// Chebyshev base count
			// FIXME some intelligent way of computing Ncb
			int psi_nodes = psi_a.getN() - psi_a.getL() - 1;
			int sturm_nodes = S_u.getN() - S_u.getL() - 1;
			int Ncb = 20 * (psi_nodes + sturm_nodes + 1);
			
			// get far radius
			Far = std::max(psi_a.far(1e-8), S_u.far(1e-8));
			
			// create the approximation
			Chebyshev cb (
				integrand,	// function to approximate
				Ncb,		// number of terms total
				0,			// left bound of interval
				Far			// right bound of interval
			);
			
			// integrate the Chebyshev approximation
			Cb = cb.integrate();
			
			// evaluate integration boundaries
			Cb_0   = Cb.approx(0,   1e-8);
			Cb_inf = Cb.approx(Far, 1e-8);
		}
		
		/**
		 * Evaluate the Phi function:
		 *     Return either 'integrate_low(x)' or 'integrate_inf(x)'
		 *     depending on the user's choice during object creation.
		 */
		double operator() (double x) const
		{
			return (this->*Eval)(x);
		}
		
		/**
		 * Evaluate proper Phi integral.
		 */
		double integrate_low (double x) const
		{
			// FIXME use Clenshaw recurrence Cb.clenshaw(...)
			
			if (x < 0)
				return 0;
			
			if (x > Far)
				return Cb_inf;
			
			return Cb.approx(x, 1e-8) - Cb_0;
		}
		
		/**
		 * Evaluate improper Phi integral.
		 */
		double integrate_inf (double x) const
		{
			// FIXME use Clenshaw recurrence Cb.clenshaw(...)
			
			if (x < 0)
				return Cb_inf;
			
			if (x > Far)
				return 0;
			
			return Cb_inf - Cb.approx(x, 1e-8);
		}
		
	private:
		
		double Far;
		Chebyshev Cb;
		double Cb_0;
		double Cb_inf;
		
		double (PhiFunction::* Eval) (double) const;
	};
	
	/**
	 * Integrate
	 * \f[
	 *     \Xi = \int_0^\inf r_1 P_\alpha(r_1) S_u^\ast(r_1) \mathrm{d}r_1
	 * \f]
	 */
	static double XiIntegral(HydrogenFunction const & psi, SturmianFunction const & S);
	
	/**
	 * Compute the integral
	 * \f[
	 *     I_{ij}^{\alpha} = \int_0^\infty \! \int_0^\infty P_{N_\alpha L_\alpha}(r_1)
	 *     \chi_{l_\alpha}(k_\alpha r_2) \left( \frac{r_<^{\lambda_\alpha}}{r_>^{\lambda_\alpha+1}}
	 *     - \frac{\delta_{\lambda_\alpha}^0}{r_2} - \delta_{\lambda_\alpha}^0
	 *     U_\alpha(r_2) \right) S_{N_i L_i}(r_1) S_{N_j L_j}(r_2) r_1 r_2 \mathrm{d}r_1
	 *      \mathrm{d}r_2 \ ,
	 * \f]
	 * which is for \f$ \lambda_\alpha = 0 \f$
	 * \f[
	 *     I_{ij}^{(\alpha)} = \int_0^\infty P_{N_\alpha L_\alpha}(r_1) S_{N_i L_i}(r_1)
	 *     \left[\int_0^{r_1}(r_2-r_1) \chi_{l_\alpha}(k_\alpha r_2) S_{N_j L_j}(r_2)
	 *     \mathrm{d}r_2 - \int_0^\infty r_1 r_2 U_\alpha(r_2) \chi_{l_\alpha}(k_\alpha
	 *     r_2) S_{N_j L_j}(r_2) \mathrm{d}r_2\right] \mathrm{d}r_1
	 * \f]
	 * and for \f$ \lambda_\alpha > 0 \f$
	 * \f[
	 *     I_{ij}^{(\alpha)} = \int_0^\infty P_{N_\alpha L_\alpha}(r_1) S_{N_i L_i}(r_1)
	 *     \left[r_1^{-\lambda_\alpha}\int_0^{r_1} r_2^{\lambda_\alpha+1} \chi_{l_\alpha}
	 *     (k_\alpha r_2) S_{N_j L_j}(r_2) \mathrm{d}r_2 + r_1^{\lambda_\alpha+1}
	 *     \int_{r_1}^\infty r_2^{-\lambda_\alpha} \chi_{l_\alpha}(k_\alpha r_2)
	 *     S_{N_j L_j}(r_2) \mathrm{d}r_2 \right] \mathrm{d}r_1 \ .
	 * \f]
	 */
	static double computeI (
		int lam_a,
		HydrogenFunction const & psi_a,
		DistortedWave const & chi_a,
		DistortingPotential const & U_a,
		SturmianFunction const & sturm1,
		SturmianFunction const & sturm2
	);
	
	/**
	 * Compute the matrix \f$ E^{(+)} -\hat{H}_{\mathrm{full}} \f$ in the Sturmian basis
	 * \f[
	 *     \left<\mathbf{r}\right| \! \left. S_i \right> =
	 *     S_{N_i L_i}(r) Y_{L_i M_i}(\mathbf{\hat{r}}) \ .
	 * \f]
	 * The implemented expression is
	 * \f[
	 *     A_{E,ijkl}^{(+)} = E^{(+)} \Sigma_{ik} \Sigma_{jl}
	 *                      - \left<\hat{H}_{\mathrm{at}}^{(1)}\right>_{ik}^{(1)} \Sigma_{jl}
	 *                      - \Sigma_{ik} \left< \hat{H}_{\mathrm{free}}^{(2)} \right>_{jl}^{(2)}
	 *                      - \left< U_{ik}^{(2)} \right>^{(2)}_{jl} \ ,
	 * \f]
	 * where
	 * \f[
	 *     \left<\hat{H}_{\mathrm{at}}^{(1)}\right>_{ik}^{(1)} = \left(
	 *     \frac{1}{2}N_k - 1\right) \delta_{ik} + \frac{1}{4}\delta_{L_i L_k}
	 *     \delta_{M_i M_k} \left[ \alpha_{N_k}^{(+)} \delta_{N_k+1,N_i}
	 *     + \alpha_{N_k}^{(-)} \delta_{N_k-1,N_i} \right] \ ,
	 * \f]
	 * \f[
	 *     \left<\hat{H}_{\mathrm{free}}^{(2)}\right>_{jl}^{(2)} = 
	 *     \frac{1}{2}N_l \delta_{jl} + \frac{1}{4}\delta_{L_j L_l}
	 *     \delta_{M_j M_l} \left[ \alpha_{N_l}^{(+)} \delta_{N_l+1,N_j}
	 *     + \alpha_{N_l}^{(-)} \delta_{N_l-1,N_j} \right] \ ,
	 * \f]
	 * \f[
	 *     U_{ik}(\mathbf{r}_2) = \int_0^\infty S_i(\mathbf{r}_1) S_k(\mathbf{r}_1)
	 *     \left( \frac{1}{r_{12}} - \frac{1}{r_2}\right) \mathrm{d}^3\!\mathbf{r}_1
	 * \f]
	 * and the overlap matrix
	 * \f[
	 *     \Sigma_{ik} = \int_0^{\infty} S_i(\mathbf{r}) S_k(\mathbf{r})
	 *     \mathrm{d}^3\!\mathbf{r} = N_k \delta_{ik} - \frac{1}{2} \delta_{L_i L_k}
	 *     \delta_{M_i M_k} \left[ \alpha_{N_k}^{(+)} \delta_{N_k+1,N_i} + \alpha_{N_k}^{(-)}
	 *     \delta_{N_k-1,N_i} \right] \ .
	 * \f]
	 * The shifting coefficients \f$ \alpha_n^{(\pm)} \f$ are (for a given \f$ \ell \f$
	 * \f[
	 *     \alpha_n^{(+)} = \sqrt{(n-\ell)(n+\ell+1)} \ ,
	 * \f]
	 * \f[
	 *     \alpha_n^{(-)} = \sqrt{(n-\ell-1)(n+\ell)} \ .
	 * \f]
	 * \warning In the present implementation, the distorting potential is zero.
	 * \param E Total energy of the system in Rydbergs.
	 */
	static double computeA (double E, int Nf1, int Nf2, int Ni1, int Ni2, int L1, int L2);
};

#endif
