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
#include <limits>
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
		) : Lam(lam), U(U), Psi(psi), Psin(psin) {}
		
		// evaluated integral
		double operator() (double x2) const
		{
			if (not finite(x2))
				return 0.;
			
			// finite integrand:
			auto integrand1 = [&](double x1) -> double {
				if (x2 == 0.)
					return 0.;
				return Psin(x1) * pow(x1/x2, Lam) * Psi(x1);
			};
			CompactIntegrand<decltype(integrand1),double> R1(integrand1, 0., Inf, 1.0);
			double i1 = (x2 == 0) ? 0. : ClenshawCurtis<decltype(R1),double>(R1, R1.scale(0.), R1.scale(x2), 1.0, 1e-10) / x2;
			
			// infinite integrand:
			auto integrand2 = [&](double x1) -> double {
				if (x1 == 0.)
					return 0.;
				return Psin(x1) * pow(x2/x1, Lam) * Psi(x1) / x1;
			};
			CompactIntegrand<decltype(integrand2),double> R2(integrand2, 0., Inf, 1.0);
			double i2 = ClenshawCurtis<decltype(R2),double>(R2, R2.scale(x2), R2.scale(Inf), 1.0, 1e-10);
			
			return i1 + i2;
		};
		
	private:
		
		/// λ (a constructor parameter)
		int Lam;
		
		/// U (a constructor parameter)
		DistortingPotential U;
		
		/// psi (a constructor parameter)
		HydrogenFunction Psi;
		
		/// psin (a constructor parameter)
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
		) : Cb_inf(0.), Lam(lam), U(U), Diag(psin == psi), 
			Zero(Diag and (DistortingPotential(psi) == U)),
			Integrand(psin,lam,U,psi),
			CompactIntegral(Integrand,0.,false,1.0)
		{
			// shorthand for infinity
			static const double inf = std::numeric_limits<double>::infinity();
			
			// the case λ = 0 is easy:
			if (lam == 0)
			{
				if (Zero)
					return;
				
				if (Diag)
				{
					auto integrand = [&](double x1) -> double {
						return (x1 == 0.) ? 0. : gsl_sf_pow_int(psi(x1), 2) / x1;
					};
					
					Cb_inf = ClenshawCurtis<decltype(integrand),double>(integrand, 0., inf, 1.0, 1e-8);
					
					return;
				}
				
				/* otherwise continue */
			}
			
			/// DEBUG
// 			std::ofstream phi;
// 			phi.open("compactphi.dat");
// 			for (int i = -100; i <= 100; i++)
// 			{
// 				double x = 0.01 * i;
// 				double ci = CompactIntegral(x);
// 				phi << x << "\t" << ci << "\n";
// 			}
// 			phi.close();
			///
			
			// convergence loop
			for (int N = 16; ; N *= 2)
			{
// 				std::cout << "[PhiFunctionDir] N = " << N << std::endl;
				
				// construct a Chebyshev approximation of the compactified function
				CompactIntegralCb = Chebyshev<double,double> (CompactIntegral, N, -1., 1.);
				
				/// DEBUG
// 				std::cout << CompactIntegralCb.str() << std::endl;
// 				std::ostringstream oss;
// 				oss << "CbN" << N << ".dat";
// 				std::ofstream ofs;
// 				ofs.open(oss.str().c_str());
// 				double cbinf = CompactIntegralCb.clenshaw(1., N-1);
// 				for (int i = -100; i <= 100; i++)
// 				{
// 					double x = 0.01 * i;
// 					ofs << x << "\t" << CompactIntegralCb.clenshaw(x,N-1) - cbinf << "\n";
// 				}
// 				ofs.close();
				///
				
				// check convergence
				if (CompactIntegralCb.tail(1e-10) != N)
					break;
			}
			
			// get optimal truncation index
			Tail = CompactIntegralCb.tail(1e-10);
			
			// evaluate Phi at positive infinity
			Cb_inf = CompactIntegralCb.clenshaw(1., Tail);
		}
		
		double operator() (double x) const
		{
			if (Lam == 0)
			{
				if (Zero) return 0.;
				if (Diag) return Cb_inf - U.plusMonopole(x);
			}
			
			return Cb_inf - CompactIntegralCb.clenshaw(CompactIntegral.scale(x), Tail);
		}
		
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

