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
 * \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_COULOMB
#define HEX_COULOMB

#include "complex.h"

//
//
// The following code has been adapted from
//    Michel N., Comp. Phys. Comm. 176 (2007) 232-249.
//
//

/**
 * \brief Direct integration of the Coulomb equation
 *
 * This class is adapted from:
 *    Michel N.: <i>Precise Coulomb wave functions for a wide range of complex ℓ, η and z</i>,
 *      Comp. Phys. Comm. 176 (2007) 232-249.
 * 
 * One uses the Burlisch-Stoer-Henrici method, where one integrates on different meshes
 * with the Henrici method, and then uses the Richardson method to get the final result by extrapolation.
 * Numerical Recipes, Chap. 16.4 .
 */
class ODE_integration
{
public:
    
    ODE_integration (Complex const & l_1, Complex const & two_eta_1);
    
    /**
     * \brief Integration of u''(r) = F(r,u(r)) with the Bulirsch-Stoer-Henrici method.
     *
     * Initials conditions : r0,u0=u(r0),du0=du/dr(r0)
     * Obtained functions : r,u=u(r),du=du/dr(r)
     *
     * See Numerical Recipes for the method.
     */
    void operator() (
        Complex const & r0,
        Complex const & u0,
        Complex const & du0,
        Complex const & r,
        Complex & u,
        Complex & du
    ) const;
    
private:
    
    /**
     * \brief Extrapolation in h=0 of a table of function values h close to h=0
     * 
     * \return Extrapolated value of the points f[h(0)]...f[h(n-1)] in h=0.
     */
    Complex extrapolation_in_zero (const unsigned int n, Complex * const f) const;
    
    /**
     * Calculation of F(z,u(z)) in u''(z) = F(z,u(z))
     * \f[
     *   F(z,u(z))=\left(\frac{l(l+1)}{z^2} + \frac{2\eta}{z} - 1\right) u(z)
     * \f]
     */
    Complex F_r_u (Complex const & r, Complex const & u) const;
    
    /**
     * \brief Integration with discretization of u''(r)=F(r,u(r)) with the Henrici method.
     * 
     * See Numerical Recipes for the method.
     *
     * Initials conditions : r0,u(r0),du/dr(r0).
     * Obtained functions : r,u(r),du/dr(r).
     */
    void integration_Henrici (
        const unsigned int m,
        Complex const & h,
        Complex const & r0,
        Complex const & u0,
        Complex const & du0,
        Complex const & r,
        Complex & u,
        Complex & du
    ) const;
    
    const Complex l,ll_plus_one;  // angular momentum,l(l+1).
    const Complex two_eta;        // 2.eta, with eta the Sommerfeld parameter.
    
    unsigned int m_tab[8];                                 // integers used in the extrapolation method.
    double one_over_m_tab[8],interpolation_term_tab[8][8]; // doubles used in the extrapolation method.
};

/**
 * \brief Class to calculate the Coulomb wave functions
 * 
 * This class is adapted from:
 *    Michel N.: <i>Precise Coulomb wave functions for a wide range of complex ℓ, η and z</i>,
 *      Comp. Phys. Comm. 176 (2007) 232-249.
 */
class Coulomb_wave_functions
{
public:
    
    Coulomb_wave_functions ();
    
    /**
     * \brief Constructor.
     *
     * Constants are defined in the constructor, 
     * plus a pointer to class ODE_integration, ODE_ptr, to integrate numerically the regular Coulomb wave function.
     * 
     * \param is_it_normalized_c true if one wants normalized functions, i.e. the standard normalization,
     *                           false if one wants F -> F/C(l,eta) and H+/H-/G -> H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1 and |z| small.
     * \param l_c orbital angular momentum.
     * \param eta_c Sommerfeld parameter.
     */
    Coulomb_wave_functions (const bool is_it_normalized_c, Complex const & l_c, Complex const & eta_c);
    
    // (shallow) copy constructor
    Coulomb_wave_functions(Coulomb_wave_functions const & WW);
    
    // destructor
    ~Coulomb_wave_functions ();
    
    /// Storage of initial conditions debut,F(debut),F'(debut)
    void F_dF_init (Complex const & z, Complex const & F, Complex const & dF);
    
    /**
     * \brief Regular wave function and derivative.
     * 
     * One calculates F(z) and F'(z), so F(z) ~ z^{l+1} for z -> 0 if is_it_normalized is false, 
     *                                   F(z) ~ C(l,eta).z^{l+1} for z -> 0 if is_it_normalized is true.
     * If |z| <= 0.5, one uses the power series.
     *
     * If |z| > 0.5 and Re[z] < 0, one calculates F(z) from F[l,-eta,-z] with F_dF_with_symmetry_relations.
     *
     * If |z| > 0.5 and Re[z] >= 0, and 1+l+/-i.eta no negative integer, one first tries the asymptotic series formula.
     * If it failed, one integrates directly the regular Coulomb wave function with F_dF_direct_integration.
     * If 1+l+/-i.eta is a negative integer, this is the only available method besides power series so it is accepted.
     * If 1+l+/-i.eta is no negative integer but it failed again, 
     * one calculates the Coulomb wave function H[omega] with H_dH_direct_integration and omega = sign(Im[z]).
     * omega is chosen so one cannot encounter the branch cut of h[omega].
     * H[-omega] is calculated from H[omega] and continued fractions h[omega] and h[-omega].
     * One then has F(z) = (H[omega] - H[-omega])/(2.i.omega.norm), F'(z) = (H'[omega] - H'[-omega])/(2.i.omega.norm),
     * with norm = 1 if the wave functions are normalized and C(l,eta)^2 if not.
     * The formula is stable as one uses this case only when |F(z)| > 0.1 .
     *
     * One takes only real parts if l, eta and z are real.
     * At the end of the function, one puts {debut,F_debut,dF_debut} equal to {z,F,dF}.
     *
     * \param z variable of the Coulomb wave function.
     * \param F regular Coulomb wave function in z and derivative in z.
     * \param dF
     */
    void F_dF (Complex const & z, Complex & F, Complex & dF) const;
    
    /**
     * \brief Calculation of G(z) and G'(z).
     * 
     * One calculates the irregular Coulomb wave function from H+ and F.
     * If 1+l+i.omega.eta is a negative integer, G is by definition H[-omega]. 
     * If not, one uses the formulas : 
     * G(z) = H+(z) - i.F(z), G'(z) = H+'(z) - i.F'(z) if is_it_normalized is true,
     * G(z) = H+(z) - i.Cl_eta^2.F(z), G'(z) = H+'(z) - i.Cl_eta^2.F'(z) if not.
     * There is no numerical inaccuracy as G is never a minimal solution.
     * One takes only real parts if l, eta and z are real.
     *
     * \param z variable of the Coulomb wave function.
     * \param G irregular Coulomb wave functions and derivatives.
     * \param dG
     */
    void G_dG (Complex const & z, Complex & G, Complex & dG) const;
    
    /**
     * \brief Calculation of H[omega](z) and H'[omega](z). 
     * 
     * One first tries the asymptotic expansion formula if 1+l+/-i.eta is no negative integer.
     * On uses logs if the unscaling factor underflows or overflows.
     * If it failed, and imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
     * one calculates H[omega](z) and H'[omega](z) with the first order expansion method.
     * If one is not in this case, one calculates F(z)and F'(z).
     * If |Im[l]| > 1 and |z| <= 1, one tries the expansion formula with F(l,eta,z) and F(-l-1,eta,z).
     * If not, or if it failed, one uses the continued fraction formula.
     * If l,eta and z are real, one rewrites H[omega] as H[omega] = Re[H[omega]] + i.omega.norm.Re[F] to avoid numerical inaccuracies for Im[H[omega]].
     * norm is 1 if is_it_normalized is true, C(l,eta)^2 if not.
     * 
     *
     * \param omega 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
     * \param z variable of the Coulomb wave function.
     * \param H : H+(z) and H+'(z) if omega=1, H-(z) and H-'(z) if omega=-1.
     * \param dH
     */
    void H_dH (const int omega, Complex const & z, Complex & H, Complex & dH) const;
    
    /**
     * \brief Calculation of the scaled H[omega](z) and H'[omega](z). 
     * 
     * They are H(omega)(z).exp[-i.omega.[z - eta.log[2z]]] and dH(omega)/dz(z).exp[-i.omega.[z - eta.log[2z]]].
     * One first tries the asymptotic expansion formula if 1+l+/-i.eta is no negative integer.
     * If it failed, and imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
     * one calculates H[omega](z) and H'[omega](z) with the first order expansion method.
     * If one is not in this case, one calculates F(z) and F'(z).
     * If |Im[l]| > 1 and |z| <= 1, one tries the expansion formula with F(l,eta,z) and F(-l-1,eta,z).
     * If not, or if it failed, one uses the continued fraction formula.
     * If l,eta and z are real, one rewrites H[omega] as H[omega] = Re[H[omega]] + I.omega.norm.Re[F], to avoid numerical inaccuracies for Im[H[omega]].
     * norm is 1 if is_it_normalized is true, C(l,eta)^2 if not.
     * One uses logs if the scaling factor underflows or overflows.
     *
     * \param omega 1 if one calculates H+(z) and H+'(z) scaled, -1 if one calculates H-(z) and H-'(z) scaled.
     * \param z variable of the Coulomb wave function.
     * \param H_scaled H[omega](z).exp[-i.omega.[z - eta.log[2z]]) and H'[omega](z).exp[-i.omega.[z - eta.log[2z]]).
     * \param dH_scaled
     */
    void H_dH_scaled (const int omega, Complex const & z, Complex & H, Complex & dH) const;
    
    const Complex l,eta; // Angular momentum and Sommerfeld parameter.
    
    const bool is_it_normalized;
    // true if F(z) ~ C(l,eta).z^{l+1} in 0, false if F(z) ~ z^{l+1} in 0.
    
private:
    
    /**
     * \brief Calculation of the asymptotic series
     * 
     * Asymptotic expansion: 
     * S(+/-)(z) = 1.0+\sum_(n=1)^N a[n] with a[n+1]=a[n].[n.[n+1+/-2i.eta]+i.eta.(i.eta+/-1)-l(l+1)]/[+/-2i.(n+1)]/z, n >= 0 and a[0] = 1.
     * This expansion diverges : it is only useful with the smallest term summation method.
     * The test of convergence is max(|a[n]|oo,|n.a[n]/z|oo), so the largest norm of the term of series of function and derivative.
     * Practically, one stops when test < precision (it worked) or when test is not finite (it failed). 
     * After that, one tests the series with the wronskian of H[omega] and H[-omega]. 
     *
     * \param omega 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
     * \param one_over_z 1/z.
     * \param sum sum[0] is the series in H(omega), sum[1] the one in H(-omega).
     * \param dsum dsum[0] is the series in H'(omega), dsum[1] the one in H'(-omega).
     * \param is_it_successful true if the series converged and |wronskian - 2.i.omega|oo is smaller than precision, false if not.
     */
    void asymptotic_series (const int omega, Complex const & one_over_z, Complex sum[], Complex dsum[], bool & is_it_successful) const;
    
    /**
     * \brief  Calculation of f = F'/F with a continued fraction.
     * 
     * One calculates the ratio f = F'/F with the continued fraction of the associated hypergeometric confluent function.
     * One uses Lentz's method.
     * One has : f = [b[0] + a[1]/b[1]+ a[2]/b[2]+ ... a[n]/b[n]+ ...]/z with :
     * b[0] = l + 1 + i.omega.z, a[n] = -2.i.omega.[1 + l + i.omega.eta] + (n-1).[-2.i.omega.z], b[n] = 2l + 2 + 2.i.omega.z + n-1.
     * omega is 1 or -1, and theoretically the result is the same.
     * If they are not equal numerically, omega = sign[-Im [z]] gives usually the best result.
     * If 1+l+i.omega.eta is a negative integer, f[omega] is finite and must be used.
     *
     * \param z Variable of the Coulomb wave function.
     * \param omega 1 or -1. Both values should be tried to test stability.
     */
    Complex continued_fraction_f (Complex const & z, const int omega) const;
    
    /**
     * \brief Calculation of h(omega) = H(omega)'/H(omega) with a continued fraction.
     *
     * One calculates the ratio h = H'/H with the continued fraction of the associated hypergeometric confluent function.
     * One uses Lentz's method.
     * One has : h = [b[0] + a[1]/b[1]+ a[2]/b[2]+ ... a[n]/b[n]+ ...].i.omega/z with :
     * b[0] = z - eta, a[n] = (1 + l + i.omega + n-1).(i.omega.eta - l + n-1), b[n] = 2[z - eta] + i.omega.n .
     *
     * If the number of iterations reaches 100000 and |z| > 0.5, the convergence is too slow. One is probably very close to the imaginary axis.
     * If l=0 and |z| <= 0.5, the direct integration is still done as it is stable.
     * If 1+l+/-i.eta is negative integer, it has to be done otherwise it is too long even if |z| is not exceedingly small.
     * One first considers Re[z] >= 0.
     * One takes the starting point z0 = 2 + i.(Im[z] + 2.sign[Im[z]]) (0.6 + i.sign[Re[z]].0.6.sign[Im[z]] if |z| <= 0.5).
     * Then, one calculates H[omega],H[omega]' at the starting point, and one integrates numerically H[omega] until z.
     * h(omega)(z) is then H[omega](z)'/H[omega](z).
     * If Re[z] < 0, one calculates H[-omega],H[-omega]' with l, eta -> -eta, and z -> -z.
     * One uses the Coulomb wave functions class cwf_minus_eta_ptr defined with l and -eta. 
     * The ratio h(omega,l,eta,z) is then equal to -h(-omega,l,-eta,-z).
     * To avoid infinite loops, continued_fraction_h must not be used in this integration. So, the starting point is chosen so |H[+/-]| is very likely
     * to increase in modulus. Then, one always integrates forward. Forward integration is enforced putting
     * is_H_dir_int_naive to true. It is put to false again at the end of the calculation.
     *
     * \param z Variable of the Coulomb wave function.
     * \param omega 1 for the outgoing wave function ratio H+'/H+, -1 for the incoming wave function ratio H-'/H-.
     */
    Complex continued_fraction_h (Complex const & z,const int omega) const;
    
    /**
     * \brief Calculation of F and F' by power series.
     * ----------------------------------------
     * It is used only when |z| <= 0.5, to avoid numerical inaccuracies.
     *
     * F(z) = norm.z^(l+1).\sum a[n], n in [0:+oo[, where :
     * a[0] = 1.0.
     * a[1] = z.eta/(l+1).
     * a[n] = (2.z.eta.a[n-1] - a[n-2].(z^2))/(n.(n+2l+1)),n >= 2.
     *
     * The z = 0 case is treated first. It is defined only for Re[l] > 0 or l = 0. The program aborts for other cases.
     * Norm is C(l,eta) if one uses normalized wave functions, 1.0 if not.
     * So, one multiplies by C(l,eta) at the end if one uses normalized functions.
     * If there is overflow or underflow for C(l,eta) in this last case, one uses logs of F,F' and C(l,eta) for the calculation.
     *
     * \param z variable of the Coulomb wave function.
     * \param F regular wave function
     * \param dF  derivative
     * The test of convergence is |(n+l-1).a[n-2]|oo + |(n+l).a[n-1]|oo, as one of the two can be zero even before convergence.
     */
    void F_dF_power_series (Complex const & z, Complex & F, Complex & dF) const;
    
    /**
     * \brief Calculation of F(z) and F'(z) with asymptotic series
     * 
     * F(z) = [H+(z) - H-(z)]/[2.i.norm]. 
     * F'(z) = [H+'(z) - H-'(z)]/[2.i.norm]. 
     * In this routine, Re[z] >= 0, so there is no branch cut problem.
     *
     * H+(z) = exp[i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].S+(z) .
     * H-(z) = exp[-i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].S-(z) .
     *
     * H+'(z) = exp[i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].[S+'(z) + S+(z).i.(1 - eta/z)] .
     * H-'(z) = exp[-i.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].[S-'(z) - S-(z).i.(1 - eta/z)] .
     *
     *
     * S+ and S- and derivatives are calculated in asymptotic_series. If is_it_successful is true, the series are meaningful. If not, one leaves the routine.
     * Norm is C(l,eta) if one uses normalized wave functions, 1.0 if not.
     * If there is overflow or underflow for C(l,eta) in this last case, one uses logs of F,F' and C(l,eta) for the calculation.
     *
     * \param z Variable of the Coulomb wave function.
     * \param F regular wave function to calculate.
     * \param dF derivative
     * \param is_it_successful true if the calculation converged, i.e. the series are good up to precision and the wronskian of H+,H- up to precision, false if not.
     */
    void asymptotic_expansion_F_dF (Complex const & z, Complex & F, Complex & dF, bool & is_it_successful) const;
    
    /**
     * \brief Calculation of H(omega) and dH(omega)/dz (scaled) with asymptotic series
     *
     * H[omega](z) = exp[i.omega.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].S[omega](z) for Re[z] >= 0.
     * H[omega]'(z) = exp[i.omega.[z - eta.log[2z] - l.Pi/2 + sigma(l,eta)]].[S[omega]'(z) + i.omega.(1 - eta/z).S[omega](z)] for Re[z] >= 0.
     *
     * S[omega](z) is the asymptotic series and S[omega]'(z) its derivative calculated in asymptotic_series.
     * If they did not converge, one leaves the routine.
     *
     * The negative cut is taken into account if log [cut_constant_AS] is finite, i.e. cut_constant_AS is not exactly zero : 
     *
     * If Re[z] < 0.0 and omega.Im[z] < 0.0 :
     * --------------------------------------
     * H[omega] = H[omega][ASd] + cut_constant_AS_plus.H[-omega][ASd].
     * H[omega][ASd] is given by directly by the asymptotic series.
     * H[-omega][ASd] is given by the asymptotic series, which gives the good result H[-omega] in this case as one is not in its bad quadrant.
     *
     * The function is scaled, so one returns : H[omega](z).exp[-i.omega.[z - eta.log[2z]]] and H[omega]'(z).exp[-i.omega.[z - eta.log[2z]]] .
     *
     * In the case of overflows or underflows, one uses logs.
     *
     * \param omega 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
     * \param one_over_z 1/z.
     * \param z Variable of the Coulomb wave function.
     * \param H_scaled H[omega](z).exp(-i.omega.[z - eta.log[2z]]) 
     * \param dH_scaled H[omega]'(z).exp(-iomega.[z - eta.log[2z]])
     * \param is_it_successful true if the asymptotic expansions converged, false it not
     */
    void asymptotic_expansion_H_dH_scaled (
        const int omega, Complex const & one_over_z,
        Complex & H_scaled, Complex & dH_scaled, bool & is_it_successful
    ) const;
    
    /**
     * \brief Calculation of F(z) and F'(z) by direct integration.
     * 
     * To calculate F(z) and F'(z), one integrates numerically F''(z) = [l(l+1)/z^2 + 2.eta/z - 1].F(z)
     * starting from debut, F_debut = F(debut) and dF_debut = F'(debut).
     * One always has Re[z] >= 0.0, so there is no branch cut problem.
     * If z = debut, the previous values are returned.
     * The starting point come from the stored values debut, F_debut and dF_debut.
     * If debut = 0, one puts debut = z/|2z| and calculates F(debut) and F'(debut) with power series.
     * The step of the integration is (z - debut)/N_num, with N_num = [|z - debut|/min (0.1,10.turning_point)] + 1.
     * The value of min (0.1,10/turning_point) gives a smaller step when turning_point increases, 
     * as calculations become there more difficult as |F| typically varies faster in this case.
     * The intermediates points are called z_aft. They go from debut to z, and (debut,F_debut,dF_debut) is put to {z_aft,F(z_aft),F'(z_aft)} at each step.
     * If |F| increases along the path, the integration is stable.
     * If it decreases, and if z_aft decreases in modulus or one does not integrate with constant argument (i.e. theta constant in z = |z|.exp[i.theta]) 
     * for Re[l] > -1, one reintegrates F(z) from debut = z/|2z|, as integration is usually stable at constant argument for Re[l] > -1.
     * Increase or decrease is known using the Taylor expansion of F near debut in z up to second order.
     * If one integrates with constant argument and |F| decreases, 
     * one integrates backwards from z_aft to debut with the knowledge of f(z_aft) = F'(z_aft)/F(z_aft).
     * F'(z_aft)/F(z_aft) is given by the continued fraction formula.
     * One then obtains by direct integration C.F(debut) and C.F'(debut).
     * Knowing F_debut, one deduces F(z) = 1/C and F'(z) = f(z_aft).F(z).
     * If 1+l+i.omega.eta is a negative integer, f[omega] is finite and is used.
     * Otherwise, f(z_aft) is calculated with omega = 1 and -1. If they are equal up to precision, f(omega) is correct and used.
     * omega is chosen so Re[-2.i.omega.z] < 0, 
     * for which the anomalous convergence phenomenon of Gautschi of f is the smallest (W. Gautschi, Math. Comp. Vol. 31 p.994).
     * If not, but |norm.F| < 0.1, one still has to use f[omega], as it is probably correct as F is the minimal solution, and also one has no other way to calculate F.
     * If |norm.F| > 0.1 in this case, one stops the procedure and F will be calculated from H+ and H-, given by direct integration and continued fraction formulae.
     * In this case, is_it_successful is put to false, and otherwise the integration worked and it is put to true.
     * Norm is 1.0 if one uses normalized wave functions, C(l,eta) if not.
     *
     * If F(z_aft) is not finite, one stops the integration and is_it_successful is put to false.
     *
     * \param z Variable of the Coulomb wave function.
     * \param F Regular wave function and derivative to calculate.
     * \param dF
     * \param is_it_successful false is the calculation is unstable, 
     *   i.e. |F| > 0.1 decreasing on the integration path and f(omega) is not equal to f(-omega) up to precision, true if not.
     */
    void F_dF_direct_integration (Complex const & z, Complex & F, Complex & dF, bool & is_it_successful) const;
    
    /**
     * \brief Calculation of H[omega](z) and H[omega]'(z) by direct integration.
     *
     *To calculate H[omega](z) and H'[omega](z), one integrates numerically H[omega]''(z) = [l(l+1)/z^2 + 2.eta/z - 1].H[omega](z)
     * starting from debut,H_debut = H[omega](debut) and dH_debut = H'[omega](debut).
     * There is no branch cut problem as Re[z] >= 0.
     * The starting point comes for the regular function from the stored values debut, F_debut and dF_debut.
     * If debut = 0, one puts debut = debut_omega = z/|2z| and calculates F(debut) and F'(debut) with power series.
     * Then, the starting point {debut_omega,H[omega](debut),H'[omega](debut)} is calculated 
     * from {debut,F(debut),F'(debut)} and the continued fraction h[omega](debut).
     * The first order expansions method is used if one is very close to the real axes of l,eta and z (Re[z] > 0).
     * The step of the integration is (z - debut)/N_num, with N_num = [|z - debut|/min (0.1,10.turning_point)] + 1.
     * The value of min (0.1,10/turning_point) gives a smaller step when turning_point increases, 
     * as calculations become there more difficult as |H[omega]| typically varies faster in this case.
     * The intermediates points are called z_aft. They go from debut to z, 
     * and (debut_omega,H_debut,dH_debut) is put to {z_aft,H[omega](z_aft),H'[omega](z_aft)} at each step.
     * If |H[omega]| increases along the path, the integration is stable.
     * If is_H_dir_int_naive is true, one has to integrate forward, as this integration is used to calculate the continued fraction.
     * If not, and if |H[omega]| decreases, 
     * one integrates backwards from z_aft to debut_omega with the knowledge of h[omega](z_aft) = H'[omega](z_aft)/H[omega](z_aft).
     * H'[omega](z_aft)/H[omega](z_aft) is given by the continued fraction formula.
     * One then obtains by direct integration C.H[omega](debut) and C.H'[omega](debut).
     * Knowing H_debut, one deduces H[omega](z) = 1/C and H'[omega](z) = f(z_aft).H[omega](z).
     * Increase or decrease is known using the Taylor expansion of H[omega] near debut_omega in z up to second order.
     *
     * If H(z_aft) is not finite, one stops the integration.
     *
     * \param omega 1 for H+,H+', -1 for H-,H-'.
     * \param z variable of the Coulomb wave function.
     * \param H wave function H[omega] to calculate.
     * \param dH derivative H'[omega] to calculate.
     */
    void H_dH_direct_integration (const int omega, Complex const & z, Complex & H, Complex & dH) const;
    
    /**
     * \brief Numerical partial derivatives according to l or eta.
     * 
     * One calculates here the partial derivatives according to l or eta
     * of the Coulomb wave functions F(x) or G(x) and of their derivatives with x F'(x) or G'(x),
     * with l_r=Re[l] and eta_r=Re[eta].
     * One considers here the parameters l_r and eta_r with the argument x.
     * For this, one uses the standard formula : df/d_chi(x) = [f(x,chi+eps) - f(x,chi-eps)]/[2.eps], 
     * with chi = l_r or eta_r and eps = prec_first_order_expansion*chi. f is either F, F', G or G'.
     * One then calculates F or G, F' or G' with chi+eps and chi-eps.
     * With them, one obtains dF/d_chi(x),dF'/d_chi(x) or dG/d_chi(x),dG'/d_chi(x), with chi = l or eta.
     * All functions are real, so double values are returned, taking the real part of complex variables.
     *
     * \param is_it_regular true if one calculates derivatives of F(z),F'(z), false if one calculates derivatives of G(z),G'(z).
     * \param is_it_eta true if one calculates the partial derivatives according to eta,
     *             false if one calculates the partial derivatives according to l.
     * \param x Re[z], z the argument of the wave function.
     * \param d_chi_Bx : partial derivatives of F and F' in l_r,eta_r,x according to chi = l or eta if is_it_regular is true, or G and G' if not.
     * \param d_chi_dBx
     */
    void partial_derivatives (const bool is_it_regular, const bool is_it_eta, const double x, double & d_chi_Bx, double & d_chi_dBx) const;
    
    /**
     * \brief Calculation of F(z),F'(z) or G(z),G'(z) with the first order expansion method.
     * 
     * When imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
     * one has to separate the calculation of the real and imaginary parts of (F,G)(z) and (F',G')(z), as they can differ by tens of orders of magnitude.
     * Re[z] < 0 is not considered, as G(z) is complex with Im[z] = Im[l] = Im[eta] = 0.
     * So, with z = x+i.y, eta = eta_r + i.eta_i and l = l_r+i.l_i, one calculates (F,G)[l_r,eta_r](x) and (F,G)'[l_r,eta_r](x).
     *
     * One considers here the parameters l_r and eta_r with the argument x, and the parameters l and eta with the argument z.
     *
     * One then has F(x),F'(x) from usual relations and one takes their real parts.
     * One also has H+(x),H+'(x) from usual relations and G(x) = Re[H+(x)],G'(x) = Re[H+'(x)].
     *
     * After that, one expands (F,G)(z) and (F',G')(z) in first order in y, eta_i and l_i and one has, up to y^2,eta_i^2 and l_i^2:
     *
     * (F,G)(z) = (F,G)(x) + i.y.d(F,G)/dx[omega](x) + i.eta_i.d(F,G)/d[eta](x) + i.l_i.d(F,G)/dl(x).
     * (F',G')[omega](z) = (F',G')[omega](x) + i.y.(F'',G'')[omega](x) + i.eta_i.d(F',G')/d[eta](x) + i.l_i.d(F',G')/dl(x).
     *
     * \param is_it_regular true if one calculates F(z),F'(z), false if one calculates G(z),G'(z).
     * \param z variable of the Coulomb wave function.
     */
    void first_order_expansions (const bool is_it_regular, Complex const & z, Complex & B, Complex & dB) const;
    
    /**
     * \brief Calculation of H[omega](z), H'[omega](z) with the first order expansion method.
     *
     * When imaginary parts of l,eta,z are much smaller than their real parts but not all zero, with Re[z] > 0,
     * one has to separate the calculation of the real and imaginary parts of H[omega] and H'[omega], as they can differ by tens of orders of magnitude.
     *
     * For that, one expands F(z),G(z),F'(z),G'(z) in first order in y, eta_i and l_i in first_order_expansions.
     *
     * Then, H[omega](z) = G(z) + i.omega.norm.F(z) and H'[omega](z) = G'(z) + i.omega.norm.F'(z),
     * with norm = 1 if the wave functions are normalized and C(l,eta)^2 if not.
     * One uses logs if norm underflows or overflows.
     *
     * \param omega 1 if one calculates H+(z),H+'(z), -1 if one calculates H-(z),H-'(z).
     * \param z variable of the Coulomb wave function.
     * \param H wave function H[omega] 
     * \param dH derivative H'[omega] to calculate.
     */
    void H_dH_from_first_order_expansions (const int omega, Complex const & z, Complex & H, Complex & dH) const;
    
    /**
     * \brief Calculation of H(omega) and H(omega)' with F(z), F'(z) and continued fractions.
     *
     * If 1+l-i.omega.eta is negative integer, one uses H[omega] = 1/F/(f - h[omega]) as it is the only available solution
     * linearly independent of F.
     * If 1+l+i.omega.eta is negative integer, H[omega] is proportional to F so the values F and F' are arbitrarily chosen for H[omega] and H'[omega].
     * 
     * If not, one calculates h[sign(Im[z])](z), and one has f = F'(z)/F(z).
     * One chooses sign(Im[z]), as h converges fastest for z in this region.
     * If |f - h[sign(Im[z])]|oo >= 1, F and H[sign] are numerically linearly independent so the continued fraction is meaningful.
     * Then, h[-sign(Im[z])](z) is not needed and is put to f, the worst value it can have. If not, h[-sign(Im[z])](z) is needed and calculated.
     * Then, h[omega] and h[-omega] take their values from h[sign(Im[z])]| and h[-sign(Im[z])]|.
     * 
     * If |f - h[omega]|oo > |f - h[-omega]|oo, one uses the continued fraction h[omega].
     * If Re[z] > 0 or sign(Im[z]) = omega (good quadrants), H[omega] = 1/F/(f - h[omega]) and H'[omega] = h[omega].H[omega].
     * If not, one has to take the branch cut into account, as one is in the bad quadrant of H[omega] : 
     * H[omega] = 1/F/(f - h[omega]) - cut_constant.F and H'[omega] = h[omega].H[omega] - cut_constant.F' .
     * cut_constant is cut_constant_CFa_plus if omega = 1, and cut_constant_CFa_minus if omega = -1.
     * If log [cut_constant] is not finite, it means that cut_constant is exactly zero if it is defined, so that there is no branch cut to consider.
     *
     * If |f - h[omega]|oo < |f - h[-omega]|oo, one uses the continued fraction h[-omega],
     * calculating H[omega] and H'[omega] using H[omega] = H[-omega] + constant.F .
     * constant is 2.i.omega.norm if Re[z] > 0 or sign(Im[z]) = -omega (good quadrants), and cut_constant if not.
     * cut_constant is cut_constant_CFb_plus if omega = 1, cut_constant_CFb_minus if omega = -1.
     * norm is 1 if is_it_normalized is true, C(l,eta)^2 if not.
     *
     * If cut_constant underflows or overflows, one uses logs of F,F' and log_cut_constant for the calculation.
     *
     * \param omega 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
     * \param z Variable of the Coulomb wave function.
    */
    void H_dH_with_F_dF_and_CF (const int omega, Complex const & z, Complex & H, Complex & dH) const;
    
    /**
     * \brief Calculation of H[omega],H'[omega] for complex l with the expansion formula.
     * 
     * When 2l is non-integer, one can expand H[omega] with F[l,eta,z] and F[-l-1,eta,z].
     * H[omega] = (exp[i.omega.chi].F - Fp)/sin (chi), H'[omega] = (exp[i.omega.chi].F' - Fp')/sin (chi) if wave functions are normalized,
     * H[omega] = (exp[i.omega.chi].F.C(l,eta)^2/sin (chi) + Fp/(2l+1), H'[omega] = (exp[i.omega.chi].F'.C(l,eta)^2/sin (chi) + Fp'/(2l+1) if not.
     * chi is sigma(l,eta) - sigma(-l-1,eta) - (l+0.5).Pi, and Fp is F(-l-1,eta).
     * Fp is calculated using a class Coulomb_wave_functions with parameters -l-1 and eta.
     * To avoid numerical imprecisions, sin (chi) is calculated with the stable formula -(2l+1).C(l,eta).C(-l-1,eta).
     * The validity of this expansion is checked with the wronskian of F and Fp, which must be correct up to precision :
     * F'.Fp - Fp'.F = sin (chi) if wave functions are normalized, Fp'.F - F'.Fp = 2l + 1 if not.
     * If the wronskians are numerically correct, one does the calculation and is_it_successful is put to true.
     * If not, one puts is_it_successful to false and quits the routine.
     * If C(l,eta)^2 underflows or overflows, one uses logs of F,F' and C(l,eta)^2 for the calculation.
     *
     * \param omega 1 if one calculates H+(z) and H+'(z), -1 if one calculates H-(z) and H-'(z).
     * \param z variable of the Coulomb wave function.
     * \param is_it_successful : false if the wronskian between F(l,eta,z) and F(-l-1,eta,z) is not equal to zero up to precision, true if not.
     */
    void H_dH_with_expansion (const int omega, Complex const & z, Complex & H, Complex & dH, bool & is_it_successful) const;
    
    /**
     * \brief Regular wave function and derivative from symmetry relations.
     *
     * If |z| > 0.5 and Re[z] < 0, one calculates F(l,eta,z),F'(l,eta,z) from F(l,-eta,-z),F'(l,-eta,-z) using the formulas :
     * F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta-i.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[-Pi.(eta-i.l)] if arg (z) > 0 and is_it_normalized is true,
     * F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta+i.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[-Pi.(eta+i.l)] if arg (z) <= 0 and is_it_normalized is true,
     * F(l,eta,z) = -F(l,-eta,-z).exp[i.Pi.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[i.Pi.l)] if arg (z) > 0 and is_it_normalized is false,
     * F(l,eta,z) = -F(l,-eta,-z).exp[-i.Pi.l)], F'(l,eta,z) = F'(l,-eta,-z).exp[-i.Pi.l)] if arg (z) <= 0 and is_it_normalized is false.
     *
     * F(l,-eta,-z) is calculated using the class cwf_minus_eta_ptr defined with (l,-eta).
     * The debut point of the class cwf_minus_eta_ptr is initialized with {debut,F_debut,dF_debut} and previous relations
     * if cwf_minus_eta_ptr->debut and -debut are different and debut non zero.
     *
     * If the normalization constant underflows or overflows, one uses logs.
     *
     * Variables
     * ---------
     * \param z Variable of the Coulomb wave function.
     * \param F Regular wave function and derivative to calculate.
     * \param dF
     */
    void F_dF_with_symmetry_relations (Complex const & z, Complex & F, Complex & dF) const;
    
    const bool neg_int_omega_one,neg_int_omega_minus_one;
    // neg_int_omega_one : true if 1+l+i.eta is negative integer, false if not.
    // neg_int_omega_minus_one : true if 1+l-i.eta is negative integer, false if not.
    
    const Complex sigma_l,log_Cl_eta,Cl_eta; // log[C(l,eta)], C(l,eta), sigma(l,eta)
    
    const Complex cut_constant_CFa_plus,cut_constant_CFa_minus,cut_constant_CFb_plus,cut_constant_CFb_minus;
    const Complex log_cut_constant_CFa_plus,log_cut_constant_CFa_minus,log_cut_constant_CFb_plus,log_cut_constant_CFb_minus;
    const Complex log_cut_constant_AS_plus,log_cut_constant_AS_minus;
    // cut constants and their logs for continued fractions (CFa and CFb) and asymptotic series (AS).
    // plus,minus is for omega = 1 or -1.
    // See functions log_cut_constant_AS_calc, log_cut_constant_CFa_calc and log_cut_constant_CFb_calc.
    
    const Complex exp_I_chi,exp_minus_I_chi,one_over_sin_chi; 
    // exp (i.chi), exp (-i.chi), 1/sin (chi) with chi = sigma(l,eta) - sigma(-l-1,eta) - (l+1/2).Pi .
    // They are used to calculate H+/H- from F(l,eta,z) and F(-l-1,eta,z) if |Im[l]| >= 1 and |z| <= 1 .
    
    const Complex sym_constant_arg_neg,sym_constant_arg_pos,log_sym_constant_arg_neg,log_sym_constant_arg_pos;
    // Multiplicative constants and their logs used in the following reflection formulas : 
    // F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta-i.l)] if arg (z) > 0 and is_it_normalized is true, so sym_constant_arg_pos = -exp[-Pi.(eta-i.l)],
    // F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta+i.l)] if arg (z) <= 0 and is_it_normalized is true, so sym_constant_arg_neg = -exp[-Pi.(eta-i.l)],
    // F(l,eta,z) = -F(l,-eta,-z).exp[i.Pi.l)] if arg (z) > 0 and is_it_normalized is false, so sym_constant_arg_pos = -exp[i.Pi.l)],
    // F(l,eta,z) = -F(l,-eta,-z).exp[-i.Pi.l)] if arg (z) <= 0 and is_it_normalized is false, so sym_constant_arg_neg = -exp[-i.Pi.l)].
    
    const double turning_point,prec_first_order_expansion; // turning_point : max (1,||eta| + sqrt[|l(l+1)| + |eta|^2]|).
    // prec_first_order_expansion : 0.1*sqrt_precision. It is the precision used for first_order_expansions.
    
    mutable bool is_H_dir_int_naive; // true if one integrates H+/H- forward without considering |H+/H-|, false if not. It is false except in continued_fraction_h.
    
    mutable Complex debut,F_debut,dF_debut;
    // Coulomb wave functions and derivative at z = debut.
    // It is used to integrate the Coulomb wave function faster, 
    // as debut is usually close to the argument of the Coulomb wave function so that the integration is quicker and more stable.
    
    mutable class ODE_integration *ODE_ptr;  // pointer to class ODE_integration to integrate numerically the Coulomb equation.
    
    mutable class Coulomb_wave_functions *cwf_real_ptr,*cwf_real_l_plus_ptr,*cwf_real_l_minus_ptr,*cwf_real_eta_plus_ptr,*cwf_real_eta_minus_ptr;
    // pointers to classes Coulomb_wave_functions of parameters (l_r,eta_r) (one has eta_r = Re[eta], eta_i = Im[eta], l_r = Re[l] and l_i = Im[l]), 
    // (l_r +/- epsilon[l],eta_r) and (l_r,eta_r +/- epsilon[eta]).
    // They are first put to zero and allocated in the program when they are needed.
    // They are used for the first order expansion method when |l_i| << 1, |eta_i| << 1 and |Im[z]| << Re[z] with Re[z] > 0.
    
    mutable class Coulomb_wave_functions *cwf_minus_eta_ptr,*cwf_lp_ptr;
    // pointers to classes Coulomb_wave_functions of parameters (l,-eta), (lp = -l-1,eta) and (l_r +/- precision,eta).
    // They are first put to zero and allocated in the program when they are needed.
    // They are used for symmetry relations : F(l,eta,z) \propto F(l,-eta,-z) and h[omega](l,eta,z) = -h[omega](l,-eta,-z)
    // and to calculate H+/H- from F(l,eta,z) and F(lp,eta,z) if |Im[l]| >= 1 and |z| <= 1.
    
    const bool noclean; 
    // whether to not execute "delete" in the constructor; corresponds to the empty constructor
    
    const bool owner;
    // whether this class is the owner of the attribute pointers
};
    
#endif
