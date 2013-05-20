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
#include <iostream>

#include "complex.h"
#include "coulomb.h"

//
//
// The following code has been adapted from
//    Michel N., Comp. Phys. Comm. 176 (2007) 232-249.
//
//

#define SIGN(a) (((a) < 0) ? (-1.) : (1.))

bool finite(Complex z)
{
	return finite(z.real()) and finite(z.imag());
}

const double precision = 1E-10, sqrt_precision = 1E-5;

//
// missing complex arithmetic
//

// Infinite norm of a complex number.
inline double inf_norm (Complex const & z)
{
	return std::max(std::abs(z.real()), std::abs(z.imag()));
}

// Logarithm of Gamma[z], z anywhere in the complex plane except in the Gamma function poles.
// ------------------------------------------------------------------------------------------
// If z is not finite or is a negative integer, the program returns an error message and stops.
// The Lanczos method is used. Precision : < 2E-10 in theory, < 1E-12 in almost every case.
// The method works for Re[z] > 0.
// If Re[z] <= 0, one uses the formula Gamma[z].Gamma[1-z] = Pi/sin (Pi.z).
// log[sin(Pi.z)] is calculated with the Kolbig method (K.S. Kolbig, Comp. Phys. Comm., Vol. 4, p.221 (1972)) : 
// If z = x+iy and y >= 0, log[sin(Pi.z)] = log[sin(Pi.eps)] - i.Pi.n, with z = n + eps so 0 <= Re[eps] < 1 and n integer.
// If y > 110, log[sin(Pi.z)] = -i.Pi.z + log[0.5] + i.Pi/2 numerically so that no overflow can occur.
// If z = x+iy and y < 0, log[Gamma(z)] = [log[Gamma(z*)]]*, so that one can use the previous formula with z*.
//
// Variables:
// ----------
// z,z_p_0p5,z_p_5p5 : argument of the Gamma function, z+0.5, z+5.5 
// sqrt_2Pi,log_Pi : sqrt(2.Pi), log(Pi).
// sum : Rational function in the Lanczos method.
// log_Gamma_z : log[Gamma(z)] value.
// c : table containing the seven coefficients in the expansion used in the Lanczos method.
// eps : z = n + eps so 0 <= Re[eps] < 1 and n integer.
// log_const : log[0.5] + i.Pi/2
Complex log_Gamma (Complex const & z)
{
	if (!finite(z)) std::cout<<"z is not finite in log_Gamma."<<std::endl, abort ();

	const double x = real (z),y = imag (z);

	if ((z == rint (x)) && (x <= 0)) std::cout<<"z is negative integer in log_Gamma."<<std::endl, abort ();

	if (x > 0.0)
	{
		const Complex z_p_0p5 = z + 0.5, z_p_5p5 = z + 5.5;
		const double sqrt_2Pi = 2.5066282746310005;
		const double c[7] = {1.000000000190015,
				7.618009172947146E+1,
				-8.650532032941677E+1,
				2.401409824083091E+1,
				-1.231739572450155,
				0.1208650973866179E-2,
				-0.5395239384953000E-5};
		
		Complex sum = c[0];
		for (int i = 1 ; i < 7 ; i++) sum += c[i]/(z + Complex(i));
		sum *= sqrt_2Pi;

		const Complex log_Gamma_z = log (sum) - log (z) + z_p_0p5*log (z_p_5p5) - z_p_5p5;

		return log_Gamma_z;
	}
	else if (y >= 0.0)
	{
		const int n = (x < rint (x)) ? (static_cast<int> (rint (x)) - 1) : (static_cast<int> (rint (x)));
		const double log_Pi = 1.1447298858494002;
		const Complex log_const(-M_LN2,M_PI_2),i_Pi(0.0,M_PI);
		const Complex eps = z - Complex(n),log_sin_Pi_z = (y > 110) ? (-i_Pi*z + log_const) : (log (sin (M_PI*eps)) - i_Pi*Complex(n));
		const Complex log_Gamma_z = log_Pi - log_sin_Pi_z - log_Gamma (1.0 - z);
		
		return log_Gamma_z;
	}
	else
		return conj (log_Gamma (conj (z)));
}

// Coulomb phase shift.
// --------------------
// It is given by the formula [Gamma[1+l+I.eta] - Gamma[1+l-I.eta]]/[2i].
// 0 is returned if 1+l+/-I.eta is a negative integer.
//
// Variables:
// ----------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta,one_over_two_I : i.eta,1/[2.i] .
// arg_plus,arg_minus : 1+l+i.eta, 1+l-i.eta.
// log_Gamma_plus,log_Gamma_minus : logs of Gamma[1+l+I.eta], Gamma[1+l-I.eta].
// sigma_l : returned result.
Complex sigma_l_calc (Complex const & l, Complex const & eta)
{
	const Complex Ieta(-imag (eta),real (eta)),one_over_two_I(0,-0.5);
	const Complex arg_plus = 1.0 + l + Ieta, arg_minus = 1.0 + l - Ieta;

	if ((rint (real (arg_plus)) == arg_plus) && (rint (real (arg_plus)) <= 0.0)) return 0.0;
	if ((rint (real (arg_minus)) == arg_minus) && (rint (real (arg_minus)) <= 0.0)) return 0.0;

	const Complex log_Gamma_plus = log_Gamma (arg_plus),log_Gamma_minus = log_Gamma (arg_minus);  
	const Complex sigma_l = (log_Gamma_plus - log_Gamma_minus)*one_over_two_I;

	return sigma_l;
}

// log of C(l,eta)
// ---------------
// It is given by the formula l*log[2] - eta.Pi/2 + (log[Gamma[1+l+I.eta]] + log[Gamma[1+l-I.eta]])/2.0 - log[Gamma[2l+2]].
// 0 is returned if 1+l+/-I.eta is a negative integer.
// 2l+2 should not be a negative integer : one has to use -l-1 instead of l in this case.
//
// Variables:
// ----------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta : i.eta .
// arg_plus,arg_minus : 1+l+i.eta, 1+l-i.eta.
// log_Gamma_plus,log_Gamma_minus,log_Gamma_2l_plus_2 : logs of Gamma[1+l+I.eta], Gamma[1+l-I.eta], Gamma[2l+2].
// log_Cl_eta : returned result.

Complex log_Cl_eta_calc (const Complex &l,const Complex &eta)
{
  const Complex Ieta(-imag (eta),real (eta));
  const Complex arg_plus = 1.0 + l + Ieta, arg_minus = 1.0 + l - Ieta; 

  if ((rint (real (arg_plus)) == arg_plus) && (rint (real (arg_plus)) <= 0.0)) return 0.0;
  if ((rint (real (arg_minus)) == arg_minus) && (rint (real (arg_minus)) <= 0.0)) return 0.0;

  const Complex log_Gamma_plus = log_Gamma (arg_plus),log_Gamma_minus = log_Gamma (arg_minus),log_Gamma_2l_plus_2 = log_Gamma (2.0*l + 2.0);
  const Complex log_Cl_eta = l*M_LN2 - M_PI_2*eta + 0.5*(log_Gamma_plus + log_Gamma_minus) - log_Gamma_2l_plus_2;

  return log_Cl_eta;
}

// Sin (chi) calculation
// ---------------------
// If 2l is integer, 0.0 is returned as chi is zero.
// If not, one calculates sin (chi) with chi = sigma(l,eta) - sigma(-l-1,eta) - (l+1/2).Pi .
// One uses the stable formula sin (chi) = -(2l+1).C(l,eta).C(-l-1,eta).
//
// Variables
// ---------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// sin_chi : sin (chi)

Complex sin_chi_calc (const Complex &l,const Complex &eta)
{
  if (rint (real (2.0*l)) == 2.0*l) return 0.0;

  const Complex sin_chi = -(2.*l+1.)*exp (log_Cl_eta_calc (l,eta) + log_Cl_eta_calc (-l-1.,eta));
  
  return sin_chi;
}

// exp (i.omega.chi) calculation.
// ------------------------------
// One calculates exp (i.omega.chi), with chi = sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi .
// If 2l is integer, 1.0 is returned as chi is zero.
// If not, one first calculates sin (chi) with the previous routine.
// If |sin (chi)| > 0.5, chi obtained with the formula sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi is stable so exp[i.omega.chi] follows directly.
// If not, one uses exp[i.omega.chi] = cos (chi) + i.omega.sin (chi), with cos (chi) = sqrt [1 - sin (chi)*sin (chi)].sign[Re[cos (chi)]],
// with chi given by sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi .
//
// Variables
// ---------
// omega : 1 or -1
// l : orbital angular momentum l.
// eta : Sommerfeld parameter = Coulomb_constant.Z.(2.mu/hbar^2)/(2k).
// I_omega : i.omega
// sin_chi : sin (chi)
// chi : sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi . 
// cos_chi : sign[Re[cos (chi)]].sqrt[1 - [sin(chi)]^2]
// exp_I_omega_chi : exp[i.omega.chi], returned result.

Complex exp_I_omega_chi_calc (const int omega,const Complex &l,const Complex &eta)
{
  if (rint (real (2.0*l)) == 2.0*l) return 1.0;

  const Complex I_omega(0,omega),sin_chi = sin_chi_calc (l,eta);
  const Complex chi = sigma_l_calc (l,eta) - sigma_l_calc (-l-1.,eta) - (l+0.5)*M_PI;

  if (std::abs (sin_chi) > 0.5)
  {
    const Complex exp_I_omega_chi = exp (I_omega*chi);

    return exp_I_omega_chi;
  }
  else
  {
    const Complex cos_chi = SIGN (real (cos (chi)))*sqrt (1.0 - sin_chi*sin_chi),exp_I_omega_chi = cos_chi + I_omega*sin_chi;

    return exp_I_omega_chi;
  }
}

// Precise evaluation of exp[z]-1 for z complex
// --------------------------------------------
// When |Re[z]| >= 1 or |Im[z]| >= 1, one uses directly the standard exp function as it is precise.
// Otherwise, numerical cancellations can occur.
// So, one uses the always stable formula exp[z]-1 = expm1(x) - 2.exp(x).sin^2(y/2) + i.exp(x).sin(y) 
// with x = Re[z] and y = Im[z]. expm1(x) gives a precise evaluation of exp(x)-1 for x double.

Complex expm1 (const Complex &z)
{
	const double x = real (z),y = imag (z);

	if ((abs (x) >= 1.0) || (abs (y) >= 1.0)) return (exp (z) - 1.0);

	const double expm1_x = expm1 (x),exp_x = 1.0 + expm1_x,sin_y_over_two = sin (0.5*y),sin_y = sin (y);

	return Complex (expm1_x - 2.0*exp_x*sin_y_over_two*sin_y_over_two,exp_x*sin_y);
}

// Precise evaluation of log[1+z] for z complex
// --------------------------------------------
// When |Re[z]| >= 1 or |Im[z]| >= 1, one uses directly the standard log function as it is precise.
// Otherwise, numerical cancellations can occur.
// So, one uses the always stable formula log[1+z] = log1p(x) + log1p([y/(1+x)]^2)/2 + i.atan2(y,1+x)
// with x = Re[z] and y = Im[z]. log1p(x) gives a precise evaluation of log(1+x) for x double.
// atan2(x,y) gives the arc tangent of y/x so it is in ]-Pi:Pi]. 

Complex log1p (const Complex &z)
{
  const double x = real (z),y = imag (z); 

  const double xp1 = 1.0 + x,abs_x = abs (x), abs_y = abs (y);

  if ((abs_x >= 1.0) || (abs_y >= 1.0)) return log (1.0 + z);
 
  const double y_over_xp1 = y/xp1;

  return Complex (log1p (x) + 0.5*log1p (y_over_xp1*y_over_xp1),atan2 (y,xp1));
}

// Cut constant log for the asymptotic series.
// -------------------------------------------
// The asymptotic series and H[omega] behave differently near the negative real axis.
// Then, if one is in the bad quadrant of H[omega], one has to take into account the cut directly.
// One is in the bad quadrant of H[omega] if Re[z] < 0.0 and sign(Im[z]) = -omega.
// 
// H[omega] = [H[omega] from asymptotic function formula] + (1 - exp(2.i.Pi.(i.eta - l.omega))).[H[-omega] from asymptotic function formula]
//
//
// The cut constant log is then log [1 - exp(2.i.Pi.(i.eta - l.omega))].
// Its returned imaginary part is not necessarily in ]-Pi:Pi]
//
//
// Variables:
// ----------
// omega : 1 or -1.
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta : i.eta
// l_int,Ieta_int : closest integers to Re[l],Re[i.eta]
// eps :  (Ieta - Ieta_int) - (l.omega - l_int.omega)
// two_I_Pi, two_I_Pi_eps : 2.i.Pi, 2.i.Pi.eps .
// log_cut_constant : returned result.

Complex log_cut_constant_AS_calc (const int omega,const Complex &l,const Complex &eta)
{
  const Complex Ieta(-imag (eta),real (eta));
  const double l_int = rint (real (l)), Ieta_int = rint (real (Ieta));
  const Complex eps = (Ieta - Ieta_int) - Complex(omega)*(l - l_int);

  const Complex two_I_Pi(0,2.0*M_PI),two_I_Pi_eps = two_I_Pi*eps;

  if (real (two_I_Pi_eps) > -0.1)
  {
    const Complex log_cut_constant = log (expm1 (-two_I_Pi_eps)) + two_I_Pi_eps;

    return log_cut_constant;
  }
  else
  {
    const Complex log_cut_constant = log1p (-exp (two_I_Pi_eps));

    return log_cut_constant;
  }
}

// Cut constant log for continued fractions : H[omega] from H[omega, not corrected] case.
// --------------------------------------------------------------------------------------
// The continued fraction has no cut on the negative real axis, whereas H[omega] has one.
// Then, if one is in the bad quadrant of H[omega], one has to take into account the cut directly.
// One is in the bad quadrant of H[omega] if Re[z] < 0.0 and sign(Im[z]) = -omega.
// 
// H[omega] = H[omega, not corrected] - cut_constant.F .
//
// The cut constant is 2i.omega.norm.(exp (2.i.Pi.[l.omega - i.eta]) - 1), and one takes its log.
// The imaginary part of the log is not necessarily in ]-Pi:Pi].
// Norm is 1.0 for normalized wave functions, C(l,eta)^2 for unnormalized wave functions.
//
//
// Variables:
// ----------  
// is_it_normalized : true if one wants normalized functions, i.e. the standard normalization,
//                    false if one wants F -> F/C(l,eta) and H+/H-/G -> H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1 and |z| small.
// omega : 1 or -1.
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta : i.eta .
// l_int,Ieta_int : closest integers to Re[l],Re[i.eta]
// eps : (l.omega - l_int.omega) - (Ieta - Ieta_int)
// log_norm : log[C(l,eta)^2] if is_it_normalized is false, 0.0 if it is true.
// two_I_Pi, two_I_Pi_eps : 2.i.Pi, 2.i.Pi.eps .
// log_two_I_omega : log[2.i.omega] = log[2] + i.omega.Pi/2 .
// log_cut_constant : returned result.

Complex log_cut_constant_CFa_calc (const bool is_it_normalized,const int omega,const Complex &l,const Complex &eta)
{
  const Complex Ieta(-imag (eta),real (eta));
  const double l_int = rint (real (l)), Ieta_int = rint (real (Ieta));
  const Complex eps = Complex(omega)*(l - l_int) - (Ieta - Ieta_int);
  const Complex log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta_calc (l,eta)) : (0.0);
  const Complex two_I_Pi(0,2.0*M_PI),two_I_Pi_eps = two_I_Pi*eps,log_two_I_omega(M_LN2,omega*M_PI_2);

  if (real (two_I_Pi_eps) < 0.1)
  {
    const Complex log_cut_constant = log_two_I_omega + log (expm1 (two_I_Pi_eps)) + log_norm;

    return log_cut_constant;
  }
  else
  {
    const Complex log_cut_constant = log_two_I_omega + log1p (-exp (-two_I_Pi_eps)) + two_I_Pi_eps + log_norm;

    return log_cut_constant; 
  }
}

// Cut constant log for continued fractions : H[omega] from H[-omega, not corrected] case.
// ---------------------------------------------------------------------------------------
// The continued fraction has no cut on the negative real axis, whereas H[-omega] has one.
// Then, if one is in the bad quadrant of H[-omega], one has to take into account the cut directly.
// One is in the bad quadrant of H[-omega] if Re[z] < 0.0 and sign(Im[z]) = omega.
// 
// H[omega] = H[-omega, not corrected] - cut_constant.F .
//
// The cut constant is 2i.omega.norm.exp (-2.i.Pi.[l.omega + i.eta]), and one takes its log.
// The returned imaginary part of the log is not necessarily in ]-Pi:Pi].
// Norm is 1.0 for normalized wave functions, C(l,eta)^2 for unnormalized wave functions.
//
//
// Variables:
// ----------
// is_it_normalized : true if one wants normalized functions, i.e. the standard normalization,
//                    false if one wants F -> F/C(l,eta) and H+/H-/G -> H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1 and |z| small.
// omega : 1 or -1.
// l : orbital angular momentum l.
// eta : Sommerfeld parameter = Coulomb_constant.Z.(2.mu/hbar^2)/(2k).
// Ieta : i.eta .
// l_int,Ieta_int : closest integers to Re[l],Re[i.eta]
// eps : (l.omega - l_int.omega) + (Ieta - Ieta_int)
// two_I_Pi : 2.i.Pi .
// log_norm : log[C(l,eta)^2] if is_it_normalized is false, 0.0 if it is true.
// log_two_I_omega : log[2.i.omega] = log[2] + i.omega.Pi/2 .
// log_cut_constant : returned result.

Complex log_cut_constant_CFb_calc (const bool is_it_normalized,const int omega,const Complex &l,const Complex &eta)
{
  const Complex Ieta(-imag (eta),real (eta));
  const double l_int = rint (real (l)), Ieta_int = rint (real (Ieta));
  const Complex eps = Complex(omega)*(l - l_int) + (Ieta - Ieta_int);
  const Complex log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta_calc (l,eta)) : (0.0);
  const Complex two_I_Pi(0,2.0*M_PI),log_two_I_omega(M_LN2,omega*M_PI_2),log_cut_constant = log_two_I_omega - two_I_Pi*eps + log_norm;

  return log_cut_constant;
}


//
// class members
//

ODE_integration::ODE_integration (Complex const & l_1, Complex const & two_eta_1)
	: l(l_1), ll_plus_one(l_1*(l_1+1.0)), two_eta(two_eta_1)
{
	for (int n = 0 ; n < 8 ; n++)
	for (int i = 0 ; i < n ; i++)
	{
		interpolation_term_tab[n][i] = 1.0;
		for (int j = 0 ; j < n ; j++)
			if (i != j)
				interpolation_term_tab[n][i] *= (i+1.0)/(i-j);
	}
	
	for (unsigned int k = 0 ; k < 8 ; k++)
		m_tab[k] = 2*(k+1);
	
	for (unsigned int k = 0 ; k < 8 ; k++)
		one_over_m_tab[k] = 1.0/static_cast<double> (m_tab[k]);
}

Complex ODE_integration::extrapolation_in_zero (const unsigned int n, Complex * const T) const
{  
	Complex f_in_zero = 0.;

	for (unsigned int i = 0 ; i < n ; i++)
		f_in_zero += interpolation_term_tab[n][i]*T[i];

	return f_in_zero;
}

Complex ODE_integration::F_r_u (Complex const & z, Complex const & u) const
{
	if (l == 0.)
		return (two_eta/z - 1.0)*u;

	const Complex one_over_z = 1.0/z;

	return ((ll_plus_one*one_over_z + two_eta)*one_over_z - 1.0)*u;
}

void ODE_integration::integration_Henrici (
	const unsigned int m, Complex const & h,
	Complex const & r0, Complex const & u0, Complex const & du0,
	Complex const & r, Complex & u, Complex & du
) const {
	
	const Complex h_square = h*h,half_h = 0.5*h;

	Complex delta = h*(du0 + half_h*F_r_u (r0,u0));
	u = u0 + delta;

	for (unsigned int i = 1 ; i < m ; i++)
	{
		delta += h_square*F_r_u (r0 + Complex(i)*h,u);
		u += delta;
	}
  
	du = delta/h + half_h*F_r_u (r,u);
}

void ODE_integration::operator() (
	Complex const & r0, Complex const & u0, Complex const & du0,
	Complex const & r, Complex & u, Complex & du
) const {
	
	if (r == r0) {u = u0; du = du0; return;}

	Complex r_debut = r0, u_debut = u0, du_debut = du0, H = r-r0, u_end[8],du_end[8],u_extrapolated_next,du_extrapolated_next;
	double test = 1.0;

	while (test > precision)
	{
		Complex H_over_m_tab[8];
		for (unsigned int k = 0 ; k < 8 ; k++) H_over_m_tab[k] = H*one_over_m_tab[k];
		const double inf_norm_half_H = inf_norm (H_over_m_tab[0]);
		
		while (inf_norm (r_debut - r) > inf_norm_half_H)
		{
			const Complex r_debut_plus_H = r_debut + H, r_end = (inf_norm (r - r_debut_plus_H) > inf_norm_half_H) ? (r_debut_plus_H) : (r);
			
			integration_Henrici (2,H_over_m_tab[0],r_debut,u_debut,du_debut,r_end,u_end[0],du_end[0]);
			integration_Henrici (4,H_over_m_tab[1],r_debut,u_debut,du_debut,r_end,u_end[1],du_end[1]);
			Complex u_extrapolated = extrapolation_in_zero (2,u_end); 
			
			unsigned int k = 2; 
			do
			{
				integration_Henrici (m_tab[k],H_over_m_tab[k],r_debut,u_debut,du_debut,r_end,u_end[k],du_end[k]);
				u_extrapolated_next = extrapolation_in_zero (++k,u_end);
				test = inf_norm (u_extrapolated/u_extrapolated_next - 1.0);
				u_extrapolated = u_extrapolated_next;
			}
			while ((test > precision) && (k < 7));
			
			r_debut += H;
			u_debut = u_extrapolated_next;
			du_debut = du_extrapolated_next = extrapolation_in_zero (k,du_end);
		}
		
		H *= 0.5;
		r_debut = r0;
		u_debut = u0;
		du_debut = du0;
	}

	u = u_extrapolated_next;
	du = du_extrapolated_next;
}

Coulomb_wave_functions::Coulomb_wave_functions()
	: is_it_normalized(false), 
	  neg_int_omega_one(false), 
	  neg_int_omega_minus_one(false), 
	  turning_point(false),
	  prec_first_order_expansion(0.),
	  noclean(true)
{
	
}

Coulomb_wave_functions::Coulomb_wave_functions (bool is_it_normalized_c, Complex const & l_c, Complex const & eta_c)
    : l (l_c),
      eta (eta_c), 
      is_it_normalized (is_it_normalized_c),
      neg_int_omega_one ((rint (real (l_c + Complex (-imag (eta_c),real (eta_c)))) == l_c + Complex (-imag (eta_c),real (eta_c))) && 
			 (rint (real (1. + l_c + Complex (-imag (eta_c),real (eta_c)))) <= 0.0)),
      neg_int_omega_minus_one ((rint (real (l_c - Complex (-imag (eta_c),real (eta_c)))) == l_c - Complex (-imag (eta_c),real (eta_c))) && 
			       (rint (real (1. + l_c - Complex (-imag (eta_c),real (eta_c)))) <= 0.0)),
      sigma_l (sigma_l_calc (l_c,eta_c)),
      log_Cl_eta (log_Cl_eta_calc (l_c,eta_c)),
      Cl_eta (exp (log_Cl_eta_calc (l_c,eta_c))),
      exp_I_chi (exp_I_omega_chi_calc (1,l_c,eta_c)),
      exp_minus_I_chi (exp_I_omega_chi_calc (-1,l_c,eta_c)),
      one_over_sin_chi (1.0/sin_chi_calc (l_c,eta_c)),
      log_cut_constant_CFa_plus (log_cut_constant_CFa_calc (is_it_normalized_c,1,l_c,eta_c)),
      log_cut_constant_CFa_minus (log_cut_constant_CFa_calc (is_it_normalized_c,-1,l_c,eta_c)),
      cut_constant_CFa_plus (exp (log_cut_constant_CFa_calc (is_it_normalized_c,1,l_c,eta_c))),
      cut_constant_CFa_minus (exp (log_cut_constant_CFa_calc (is_it_normalized_c,-1,l_c,eta_c))),
      log_cut_constant_CFb_plus (log_cut_constant_CFb_calc (is_it_normalized_c,1,l_c,eta_c)),
      log_cut_constant_CFb_minus (log_cut_constant_CFb_calc (is_it_normalized_c,-1,l_c,eta_c)),
      log_cut_constant_AS_plus (log_cut_constant_AS_calc (1,l_c,eta_c)),
      log_cut_constant_AS_minus (log_cut_constant_AS_calc (-1,l_c,eta_c)),
      cut_constant_CFb_plus (exp (log_cut_constant_CFb_calc (is_it_normalized_c,1,l_c,eta_c))),
      cut_constant_CFb_minus (exp (log_cut_constant_CFb_calc (is_it_normalized_c,-1,l_c,eta_c))),
      log_sym_constant_arg_neg ((is_it_normalized_c) ? (-M_PI*(eta_c+(l_c+1.)*Complex (0.0,1.0))) : (-M_PI*(l_c+1.)*Complex (0.0,1.0))),
      log_sym_constant_arg_pos ((is_it_normalized_c) ? (-M_PI*(eta_c-(l_c+1.)*Complex (0.0,1.0))) : (M_PI*(l_c+1.)*Complex (0.0,1.0))),
      sym_constant_arg_neg ((is_it_normalized_c) ? (exp (-M_PI*(eta_c+(l_c+1.)*Complex (0.0,1.0)))) : (exp (-M_PI*(l_c+1.)*Complex (0.0,1.0)))),
      sym_constant_arg_pos ((is_it_normalized_c) ? (exp (-M_PI*(eta_c-(l_c+1.)*Complex (0.0,1.0)))) : (exp (M_PI*(l_c+1.)*Complex (0.0,1.0)))), 
      turning_point (std::max (1.0,std::abs (eta_c) + sqrt (std::abs (l_c*(l_c+1.0)) + std::abs (eta_c*eta_c)))),
      is_H_dir_int_naive (false),cwf_real_ptr (0),cwf_real_eta_plus_ptr (0),cwf_real_eta_minus_ptr (0),cwf_real_l_plus_ptr (0),cwf_real_l_minus_ptr (0),
      cwf_minus_eta_ptr (0),cwf_lp_ptr (0),prec_first_order_expansion (0.1*sqrt_precision),noclean(false)
{
	ODE_ptr = new class ODE_integration (l,2.0*eta);
	debut = 0.0;
	if (real (l) >= 0.0)
	{
		F_debut = 0.0;
		dF_debut = (l == 0.) ? ((is_it_normalized) ? (Cl_eta) : (1.0)) : (0.0);
	}
}

Coulomb_wave_functions::~Coulomb_wave_functions()
{
	if (noclean)
		return;
	
	delete cwf_real_ptr;

	delete cwf_real_l_plus_ptr;
	delete cwf_real_l_minus_ptr;

	delete cwf_real_eta_plus_ptr;
	delete cwf_real_eta_minus_ptr;

	delete cwf_minus_eta_ptr;
	delete cwf_lp_ptr;

	delete ODE_ptr;
}

Complex Coulomb_wave_functions::continued_fraction_h (Complex const & z, const int omega) const
{ 
	const double small = 1E-50,large = 1E50,abs_z = std::abs (z);
	const Complex I_omega(0.0,omega),two_I_omega(0.0,2.0*double(omega)),I_omega_eta = I_omega*eta;
	const Complex a = I_omega_eta + l + 1.0,c = I_omega_eta - l,z_minus_eta = z - eta,two_z_minus_eta = 2.0*z_minus_eta;

	Complex b0 = z_minus_eta,hn = (b0 != 0.0) ? (b0) : (1E-50), Cn = hn, Dn = 0.0;
	int n = 1;
	double test;
	do
	{
		const int nm1 = n-1;
		const Complex an = (a + Complex(nm1))*(c + Complex(nm1)),bn = two_z_minus_eta + double(n)*two_I_omega,bn_plus_an_Dn = bn + an*Dn,bn_plus_an_over_Cn = bn + an/Cn;
		
		Dn = (bn_plus_an_Dn != 0.0) ? (1.0/bn_plus_an_Dn) : (large);
		Cn = (bn_plus_an_over_Cn != 0.0) ? (bn_plus_an_over_Cn) : (small);
		
		const Complex Delta_n = Dn*Cn;
		hn *= Delta_n;
		test = inf_norm (1.0 - Delta_n);
		
		if ((n++ > 100000) && ((l == 0.0) || (abs_z > 0.5) || neg_int_omega_one || neg_int_omega_minus_one))
		{
			if ((real (z) < 0.0) && (cwf_minus_eta_ptr == 0)) cwf_minus_eta_ptr = new class Coulomb_wave_functions (is_it_normalized,l,-eta);
			class Coulomb_wave_functions const & cwf = (real (z) < 0.0) ? (*cwf_minus_eta_ptr) : (*this);
			
			const Complex z00(2.0, SIGN(real(z))*(imag (z) + 2.0*SIGN (imag (z)))),z01(0.6,0.6*SIGN (real (z))*SIGN (imag (z)));
			const Complex z_start = (abs_z > 0.5) ? (z00) : (z01),debut_cwf = cwf.debut,F_debut_cwf = cwf.F_debut,dF_debut_cwf = cwf.dF_debut;
			Complex F_start,dF_start,H,dH;
			cwf.F_dF (z_start,F_start,dF_start);
			cwf.is_H_dir_int_naive = true, cwf.H_dH_direct_integration (SIGN (real (z))*omega,SIGN (real (z))*z,H,dH), cwf.is_H_dir_int_naive = false;
			cwf.debut = debut_cwf, cwf.F_debut = F_debut_cwf, cwf.dF_debut = dF_debut_cwf;
			return (SIGN (real (z))*dH/H);
		}
	}
	while (test > 1E-15);

	const Complex h = hn*I_omega/z;
	return h;
}

void Coulomb_wave_functions::asymptotic_expansion_H_dH_scaled (const int omega,Complex const & one_over_z,
							       Complex &H_scaled,Complex &dH_scaled,bool &is_it_successful) const
{  
	Complex sum[2],dsum[2];

	asymptotic_series (omega,one_over_z,sum,dsum,is_it_successful);
	if (!is_it_successful) return;

	const Complex I_omega(0,omega),two_I_omega(0,2*omega),I_omega_one_minus_eta_over_z = I_omega*(1.0 - eta*one_over_z);

	const Complex phase_shift = I_omega*(sigma_l - l*M_PI_2),exp_phase_shift = exp (phase_shift);

	H_scaled = sum[0]*exp_phase_shift;
	dH_scaled = (dsum[0] + sum[0]*I_omega_one_minus_eta_over_z)*exp_phase_shift;

	const Complex log_cut_constant_AS = (omega == 1) ? (log_cut_constant_AS_plus) : (log_cut_constant_AS_minus);

	if (one_over_z != 0.0)
	{
		const Complex z = 1.0/one_over_z;
	
		if (finite (log_cut_constant_AS) && (real (z) < 0.0) && (SIGN (imag (z)) == -omega))
		{
			const Complex factor = exp (-two_I_omega*(z - eta*(M_LN2 + log (z))) - phase_shift + log_cut_constant_AS);

			H_scaled += sum[1]*factor;
			dH_scaled += (dsum[1] - sum[1]*I_omega_one_minus_eta_over_z)*factor;
		}
	}

	if (!is_it_normalized)
	{
		if ((Cl_eta == 0.0) || (!finite (Cl_eta))) 
			H_scaled = exp (log_Cl_eta + log (H_scaled)),dH_scaled = exp (log_Cl_eta + log (dH_scaled));
		else 
			H_scaled *= Cl_eta,dH_scaled *= Cl_eta;
	}
}

void Coulomb_wave_functions::H_dH_direct_integration (const int omega,Complex const & z,Complex &H,Complex &dH) const
{ 
	const double x = real (z),y = imag (z),/*l_r = real (l),*/l_i = imag (l),/*eta_r = real (eta),*/eta_i = imag (eta);

	if (debut == 0.0) 
	{
		debut = 0.5*z/std::abs (z);
		F_dF_power_series (debut,F_debut,dF_debut);
	}

	Complex debut_omega = debut,H_debut,dH_debut;

	if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
		&& (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
		&& (!neg_int_omega_one && !neg_int_omega_minus_one))
		first_order_expansions (omega,debut_omega,H_debut,dH_debut);
	else
		H_dH_with_F_dF_and_CF (omega,debut_omega,H_debut,dH_debut);

	const class ODE_integration &ODE = *ODE_ptr;
	const double step_abs = std::min(0.1,10.0/turning_point);
	const unsigned int N_num = static_cast<unsigned int> ((std::abs (z-debut_omega)/step_abs) + 1);
	const Complex step_num = (z - debut_omega)/static_cast<double> (N_num),ll_plus_one = l*(l+1.0),two_eta = 2.0*eta;

	for (unsigned int i = N_num-1 ; i <= N_num ; i--)
	{
		const Complex z_aft = z - double(i)*step_num,one_over_debut = 1.0/debut,log_H_debut_der = dH_debut/H_debut;
		const Complex d2H_debut_over_H_debut = (ll_plus_one*one_over_debut + two_eta)*one_over_debut - 1.0;
	
		if (is_H_dir_int_naive)
			ODE (debut_omega,H_debut,dH_debut,z_aft,H,dH);
		else if (std::abs (1.0 + step_num*(log_H_debut_der + 0.5*step_num*d2H_debut_over_H_debut)) < 1.0)
		{
			const Complex h = continued_fraction_h (z_aft,omega);      
			Complex H_debut_not_normed,dH_debut_not_normed;

			ODE (z_aft,1.0,h,debut_omega,H_debut_not_normed,dH_debut_not_normed);
			H = H_debut/H_debut_not_normed;
			dH = h*H; 
		}
		else ODE (debut_omega,H_debut,dH_debut,z_aft,H,dH);

		if (!finite (H) || !finite (dH))
			std::cout<<"Numerical failure encountered in H_dH_direct_integration."<<std::endl,exit (1);

		debut_omega = z_aft,H_debut = H,dH_debut = dH;
	}
}

void Coulomb_wave_functions::H_dH_from_first_order_expansions (const int omega,Complex const & z,Complex &H,Complex &dH) const
{
	Complex F,dF,G,dG;
	first_order_expansions (true,z,F,dF);
	first_order_expansions (false,z,G,dG);

	const Complex I_omega(0,omega),norm_functions = (!is_it_normalized) ? (Cl_eta*Cl_eta) : (1.0);

	if ((norm_functions == 0.0) || (!finite (norm_functions))) 
	{
		const Complex log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta) : (0.0);

		H = G + I_omega*exp (log (F) + log_norm);
		dH = dG + I_omega*exp (log (dF) + log_norm);
	}
	else
	{
		H = G + I_omega*norm_functions*F;
		dH = dG + I_omega*norm_functions*dF;
	}
}

void Coulomb_wave_functions::H_dH_with_F_dF_and_CF (const int omega,Complex const & z,Complex &H,Complex &dH) const
{  
	Complex F,dF;
	F_dF (z,F,dF);

	const double x = real (z),y = imag (z);
	const Complex f = dF/F,two_I_omega(0.0,2.0*omega);

	if (((neg_int_omega_one && (omega == -1))) || ((neg_int_omega_minus_one && (omega == 1))))
	{
		const Complex h_omega = continued_fraction_h (z,omega);
		H = 1.0/(F*(f - h_omega)),dH =  H*h_omega;
	}
	else if (neg_int_omega_one || neg_int_omega_minus_one)
		H = F,dH = dF;
	else
	{
		const Complex h_sign = continued_fraction_h (z,SIGN (y)),h_minus_sign = (std::abs (f - h_sign) < 1.0) ? (continued_fraction_h (z,-SIGN (y))) : (f);
		const Complex h_omega = (omega == SIGN (y)) ? (h_sign) : (h_minus_sign),h_minus_omega = (omega == SIGN (y)) ? (h_minus_sign) : (h_sign);

		if (inf_norm (f - h_omega) > inf_norm (f - h_minus_omega))
		{
			H = 1.0/(F*(f - h_omega)),dH =  H*h_omega;

			const Complex log_cut_constant_CFa = (omega == 1) ? (log_cut_constant_CFa_plus) : (log_cut_constant_CFa_minus);

			if ((finite (log_cut_constant_CFa)) && (x < 0.0) && (SIGN (y) == -omega))
			{
				const Complex cut_constant_CFa = (omega == 1) ? (cut_constant_CFa_plus) : (cut_constant_CFa_minus); 
				
				if ((cut_constant_CFa == 0.0) || (!finite (cut_constant_CFa)))
					H -= exp (log_cut_constant_CFa + log (F)),dH -= exp (log_cut_constant_CFa + log (dF));
				else 
					H -= cut_constant_CFa*F,dH -= cut_constant_CFa*dF;
			}
		}
		else
		{
			const Complex H_minus_omega = 1.0/(F*(f - h_minus_omega)),dH_minus_omega = H_minus_omega*h_minus_omega;
			const Complex norm_functions = (!is_it_normalized) ? (Cl_eta*Cl_eta) : (1.0);
			const Complex cut_constant_CFb = (omega == 1) ? (cut_constant_CFb_plus) : (cut_constant_CFb_minus); 
			const Complex constant = ((x < 0.0) && (SIGN (y) == omega)) ? (cut_constant_CFb) : (two_I_omega*norm_functions);

			if ((constant == 0.0) || (!finite (constant))) 
			{
				const Complex log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta) : (0.0),log_two_I_omega(M_LN2,omega*M_PI_2);
				const Complex log_cut_constant_CFb = (omega == 1) ? (log_cut_constant_CFb_plus) : (log_cut_constant_CFb_minus);
				const Complex log_constant = ((x < 0.0) && (SIGN (y) == omega)) ? (log_cut_constant_CFb) : (log_two_I_omega + log_norm);
				
				H = exp (log_constant + log (F)) + H_minus_omega,dH = exp (log_constant + log (dF)) + dH_minus_omega;
			}
			else H = constant*F + H_minus_omega,dH = constant*dF + dH_minus_omega;
		}
	}
}

void Coulomb_wave_functions::H_dH_with_expansion (const int omega,Complex const & z,Complex &H,Complex &dH,bool &is_it_successful) const
{
	if (cwf_lp_ptr == 0) cwf_lp_ptr = new class Coulomb_wave_functions (is_it_normalized,-l-1.,eta);

	Complex F,dF,Fp,dFp;
	F_dF (z,F,dF);
	cwf_lp_ptr->F_dF (z,Fp,dFp);
		
	const Complex exp_I_omega_chi = (omega == 1) ? (exp_I_chi) : (exp_minus_I_chi);
		
	if (is_it_normalized)
	{ 
		if (inf_norm ((dFp*F - dF*Fp)*one_over_sin_chi - 1.0) < precision)
		{
			H = (exp_I_omega_chi*F - Fp)*one_over_sin_chi;
			dH = (exp_I_omega_chi*dF - dFp)*one_over_sin_chi;
		}
		else {is_it_successful = false; return;}
	}
	else  
	{
		const Complex one_over_2lp1 = 1.0/(2.*l+1.);

		if (inf_norm ((dF*Fp - dFp*F)*one_over_2lp1 - 1.0) < precision)
		{
			const Complex Cl_eta_2 = Cl_eta*Cl_eta,exp_I_omega_chi_over_sin_chi = exp_I_omega_chi*one_over_sin_chi;
			const Complex F_Cl_eta_2 = ((Cl_eta_2 == 0.0) || (!finite (Cl_eta_2))) ? (exp (2.0*log_Cl_eta + log (F))) : (Cl_eta_2*F);
			const Complex dF_Cl_eta_2 = ((Cl_eta_2 == 0.0) || (!finite (Cl_eta_2))) ? (exp (2.0*log_Cl_eta + log (dF))) : (Cl_eta_2*dF);
			
			H = exp_I_omega_chi_over_sin_chi*F_Cl_eta_2 + Fp*one_over_2lp1;
			dH = exp_I_omega_chi_over_sin_chi*dF_Cl_eta_2 + dFp*one_over_2lp1;
		}
		else {is_it_successful = false; return;}
	}
	
	is_it_successful = true;
}

void Coulomb_wave_functions::F_dF_power_series (Complex const & z,Complex &F,Complex &dF) const
{
	if (z == 0.0)
	{
		if (l == 0.) F = 0.0,dF = (is_it_normalized) ? (Cl_eta) : (1.0);
		else if (real (l) > 0) F = dF = 0.0;
		else std::cout << "F(z=0) and/or F'(z=0) are undefined." << std::endl, abort ();
	}
	else
	{
		const Complex z_square = z*z,z_two_eta = 2.0*eta*z;

		int n = 2;
		Complex an_minus_two = 1.0,an_minus_one = z*eta/(l+1.0);
	
		F = an_minus_two + an_minus_one;
		dF = (l+1.0)*an_minus_two + (l+2.0)*an_minus_one;

		while (inf_norm (an_minus_two*(double(n)+l-1.0)) + inf_norm (an_minus_one*(double(n)+l)) > precision)
		{
			const Complex an = (z_two_eta*an_minus_one - an_minus_two*z_square)/(double(n)*(double(n) + l + l + 1.0));

			F += an;
			dF += an*(double(n) + l + 1.0);
			
			n++;
			an_minus_two = an_minus_one;
			an_minus_one = an;
		}
		
		const Complex z_pow_l_plus_one = pow (z,l+1.0);
		F *= z_pow_l_plus_one;
		dF *= z_pow_l_plus_one/z; 
		
		if (is_it_normalized)
		{
			if ((Cl_eta == 0.0) || (!finite (Cl_eta))) F = exp (log_Cl_eta + log (F)),dF = exp (log_Cl_eta + log (dF));
			else F *= Cl_eta,dF *= Cl_eta;
		}
	}
}

Complex Coulomb_wave_functions::continued_fraction_f (Complex const & z,const int omega) const
{
	const double small = 1E-50,large = 1E50;

	const Complex I_omega(0.0,omega);
	const Complex a = I_omega*eta + l + 1.0,b = 2.*l + 2.;
	const Complex minus_two_I_omega_z = -2.0*I_omega*z,minus_two_I_omega_a_z = minus_two_I_omega_z*a,b_plus_two_I_omega_z = b - minus_two_I_omega_z;

	const Complex b0 = l + 1.0 + I_omega*z;
	Complex fn = (b0 != 0.0) ? (b0) : (small), Cn = fn, Dn = 0.0;  

	int n = 1;
	double test;
	do
	{
		const int nm1 = n-1;
		const Complex an = minus_two_I_omega_a_z + double(nm1)*minus_two_I_omega_z;
		const Complex bn = b_plus_two_I_omega_z + double(nm1);
		
		const Complex bn_plus_an_Dn = bn + an*Dn,bn_plus_an_over_Cn = bn + an/Cn;

		Dn = (bn_plus_an_Dn != 0.0) ? (1.0/bn_plus_an_Dn) : (large);
		Cn = (bn_plus_an_over_Cn != 0.0) ? (bn_plus_an_over_Cn) : (small);

		const Complex Delta_n = Dn*Cn;
		fn *= Delta_n;
		test = inf_norm (1.0 - Delta_n);
		n++;
	}
	while (test > 1E-15);

	const Complex f = fn/z;

	return f;
}

void Coulomb_wave_functions::asymptotic_expansion_F_dF (Complex const & z,Complex &F,Complex &dF,bool &is_it_successful) const
{
	Complex sum[2],dsum[2];

	const Complex one_over_z = 1.0/z;

	asymptotic_series (1,one_over_z,sum,dsum,is_it_successful);
	if (!is_it_successful) return;

	const Complex I(0,1),one_over_two_I(0.0,-0.5),I_one_minus_eta_over_z = I*(1.0 - eta*one_over_z);

	const Complex exp_phase_shift_plus = exp (I*(z - eta*(M_LN2 + log (z)) - l*M_PI_2 + sigma_l));
	const Complex exp_phase_shift_minus = 1.0/exp_phase_shift_plus;

	const Complex H_plus = sum[0]*exp_phase_shift_plus,dH_plus = (dsum[0] + sum[0]*I_one_minus_eta_over_z)*exp_phase_shift_plus;
	const Complex H_minus = sum[1]*exp_phase_shift_minus,dH_minus = (dsum[1] - sum[1]*I_one_minus_eta_over_z)*exp_phase_shift_minus;

	F = (H_plus - H_minus)*one_over_two_I;
	dF = (dH_plus - dH_minus)*one_over_two_I;

	if (!is_it_normalized)
	{
		if ((Cl_eta == 0.0) || (!finite (Cl_eta))) F = exp (log (F) - log_Cl_eta),dF = exp (log (dF) - log_Cl_eta);
		else F /= Cl_eta,dF /= Cl_eta;
	}
}

void Coulomb_wave_functions::F_dF_direct_integration (Complex const & z,Complex &F,Complex &dF,bool &is_it_successful) const
{
	if (z == debut) {F = F_debut,dF = dF_debut,is_it_successful = true; return;}
	if (debut == 0.0) debut = 0.5*z/std::abs (z),F_dF_power_series (debut,F_debut,dF_debut);

	const class ODE_integration &ODE = *ODE_ptr;
	const double step_abs = std::min(0.1,10.0/turning_point);
	const unsigned int N_num = static_cast<unsigned int> (rint (std::abs (z-debut)/step_abs) + 1); 
	const Complex step_num = (z - debut)/static_cast<double> (N_num),ll_plus_one = l*(l+1.0),two_eta = 2.0*eta;
	const Complex norm_functions = (!is_it_normalized) ? (Cl_eta*Cl_eta) : (1.0);

	for (unsigned int i = N_num-1 ; i <= N_num ; i--)
	{
		const Complex z_aft = z - double(i)*step_num,one_over_debut = 1.0/debut,log_F_debut_der = dF_debut/F_debut;
		const Complex d2F_debut_over_F_debut = (ll_plus_one*one_over_debut + two_eta)*one_over_debut - 1.0;

		if (std::abs (1.0 + step_num*(log_F_debut_der + 0.5*step_num*d2F_debut_over_F_debut)) < 1.0)
		{
			const Complex ratio = debut/z; 
			if ((real (l) > -1.0) && ((std::abs (imag (ratio)) > precision) || (real (ratio) > 1.0)))
			{debut = 0.0, F_dF_direct_integration (z,F,dF,is_it_successful); return;}

			Complex F_debut_not_normed,dF_debut_not_normed;

			if (neg_int_omega_one)
			{
				const Complex fp = continued_fraction_f (z_aft,1);
				ODE (z_aft,1.0,fp,debut,F_debut_not_normed,dF_debut_not_normed);
				F = F_debut/F_debut_not_normed; 
				dF = fp*F;
			}
			else if (neg_int_omega_minus_one)
			{
				const Complex fm = continued_fraction_f (z_aft,-1);
				ODE (z_aft,1.0,fm,debut,F_debut_not_normed,dF_debut_not_normed); 
				F = F_debut/F_debut_not_normed; 
				dF = fm*F;
			}
			else
			{
				const Complex f_omega = continued_fraction_f (z_aft,SIGN(-imag(z_aft))),f_minus_omega = continued_fraction_f (z_aft,-SIGN(-imag(z_aft)));

				//// Comment the following line if you want to accept the value of f(omega) = F'/F inconditionally
				if ((std::abs (F*norm_functions) > 0.1) && (std::abs (f_minus_omega/f_omega - 1.0) > precision)) {is_it_successful = false; return;}

				ODE (z_aft,1.0,f_omega,debut,F_debut_not_normed,dF_debut_not_normed);
				F = F_debut/F_debut_not_normed; 
				dF = f_omega*F;
			}
		}
		else ODE (debut,F_debut,dF_debut,z_aft,F,dF);

		debut = z_aft,F_debut = F,dF_debut = dF;

		if (!finite (F) || !finite (dF)) std::cout<<"Numerical failure encountered in F_dF_direct_integration."<<std::endl,exit (1);
	}
	is_it_successful = true;
}

void Coulomb_wave_functions::F_dF_with_symmetry_relations (Complex const & z,Complex &F,Complex &dF) const
{
	if (cwf_minus_eta_ptr == 0) cwf_minus_eta_ptr = new class Coulomb_wave_functions (is_it_normalized,l,-eta);

	if ((debut != 0.0) && (cwf_minus_eta_ptr->debut != -debut))
	{
		const double arg_debut = arg (debut);
		const Complex sym_constant_debut = (arg_debut <= 0.0) ? (1.0/sym_constant_arg_neg) : (1.0/sym_constant_arg_pos);

		cwf_minus_eta_ptr->debut = -debut; 

		if ((sym_constant_debut == 0.0) || (!finite (sym_constant_debut))) 
		{
			const Complex log_sym_constant_debut = (arg (debut) <= 0.0) ? (-log_sym_constant_arg_neg) : (-log_sym_constant_arg_pos);

			cwf_minus_eta_ptr->F_debut = exp (log_sym_constant_debut + log (F_debut)); 
			cwf_minus_eta_ptr->dF_debut = -exp (log_sym_constant_debut + log (dF_debut));
		} 
		else
		{
			cwf_minus_eta_ptr->F_debut = F_debut*sym_constant_debut; 
			cwf_minus_eta_ptr->dF_debut = -dF_debut*sym_constant_debut;
		}
	}

	const double arg_z = arg (z);
	const Complex sym_constant = (arg_z <= 0.0) ? (sym_constant_arg_neg) : (sym_constant_arg_pos);

	cwf_minus_eta_ptr->F_dF (-z,F,dF);

	if ((sym_constant == 0.0) || (!finite (sym_constant))) 
	{
		const Complex log_sym_constant = (arg_z <= 0.0) ? (log_sym_constant_arg_neg) : (log_sym_constant_arg_pos);

		F = exp (log_sym_constant + log (F)); 
		dF = -exp (log_sym_constant + log (dF));
	} 
	else 
	{
		F *= sym_constant; 
		dF *= -sym_constant;
	}
}

void Coulomb_wave_functions::asymptotic_series (const int omega,Complex const & one_over_z,
						Complex sum[],Complex dsum[],bool &is_it_successful) const
{
	sum[0] = sum[1] = 1.0;
	dsum[0] = dsum[1] = 0.0;

	double test;

	for (unsigned int i = 0 ; i <= 1 ; i++)
	{
		const double sign = (i == 0) ? (omega) : (-omega); 
		const Complex Ieta(-imag (eta),real (eta)),two_I_eta_sign = 2.0*sign*Ieta,Ieta_Ieta_plus_sign_minus_ll_plus_one = Ieta*(Ieta + sign) - l*(l+1.0);
		
		int n = 0;
		Complex an_sign = 1.0;

		do
		{
			const double n_plus_one = n + 1.0;
			const Complex sign_one_over_two_I_n_plus_one(0,-sign*0.5/n_plus_one);

			an_sign *= one_over_z*(double(n)*(n_plus_one + two_I_eta_sign) + Ieta_Ieta_plus_sign_minus_ll_plus_one)*sign_one_over_two_I_n_plus_one;

			const Complex sum_term = an_sign,dsum_term = n_plus_one*an_sign*one_over_z;

			sum[i] += sum_term;
			dsum[i] -= dsum_term;

			test = std::max (inf_norm (sum_term),inf_norm (dsum_term));
			n++;
		}
		while ((test > precision) && (finite (test)));

		if (!finite (test)) {is_it_successful = false; return;}
	}

	const Complex two_I_omega(0.0,2.0*omega);
	is_it_successful = (inf_norm (sum[1]*dsum[0] - sum[0]*dsum[1] + two_I_omega*(1.0 - eta*one_over_z)*sum[0]*sum[1] - two_I_omega) < precision);
}

void Coulomb_wave_functions::partial_derivatives (const bool is_it_regular,const bool is_it_eta,const double x,double &d_chi_Bx,double &d_chi_dBx) const
{
	const double l_r = real (l),eta_r = real (eta),chi_r = (is_it_eta) ? (eta_r) : (l_r);
	const double chi_r_plus = (chi_r != 0.0) ? (chi_r*(1.0 + prec_first_order_expansion)) : (prec_first_order_expansion);
	const double chi_r_minus = (chi_r != 0.0) ? (chi_r*(1.0 - prec_first_order_expansion)) : (-prec_first_order_expansion);

	if (is_it_eta)
	{ 
		if (cwf_real_eta_plus_ptr == 0) cwf_real_eta_plus_ptr = new class Coulomb_wave_functions (is_it_normalized,l_r,chi_r_plus);
		if (cwf_real_eta_minus_ptr == 0) cwf_real_eta_minus_ptr = new class Coulomb_wave_functions (is_it_normalized,l_r,chi_r_minus);
	}
	else
	{
		if (cwf_real_l_plus_ptr == 0) cwf_real_l_plus_ptr = new class Coulomb_wave_functions (is_it_normalized,chi_r_plus,eta_r);
		if (cwf_real_l_minus_ptr == 0) cwf_real_l_minus_ptr = new class Coulomb_wave_functions (is_it_normalized,chi_r_minus,eta_r);
	}
	
	class Coulomb_wave_functions &cwf_plus = (is_it_eta) ? (*cwf_real_eta_plus_ptr) : (*cwf_real_l_plus_ptr);
	class Coulomb_wave_functions &cwf_minus = (is_it_eta) ? (*cwf_real_eta_minus_ptr) : (*cwf_real_l_minus_ptr);

	Complex A_plus,dA_plus,A_minus,dA_minus;

	if (is_it_regular)
	{
		cwf_plus.F_dF (x,A_plus,dA_plus);
		cwf_minus.F_dF (x,A_minus,dA_minus);
	}
	else
	{
		cwf_plus.H_dH (1,x,A_plus,dA_plus);
		cwf_minus.H_dH (1,x,A_minus,dA_minus);
	}

	const double B_plus = real (A_plus),B_minus = real (A_minus),dB_plus = real (dA_plus),dB_minus = real (dA_minus);

	d_chi_Bx = (B_plus - B_minus)/(chi_r_plus - chi_r_minus);
	d_chi_dBx = (dB_plus - dB_minus)/(chi_r_plus - chi_r_minus);
}

void Coulomb_wave_functions::first_order_expansions (const bool is_it_regular,Complex const & z,Complex &B,Complex &dB) const
{
	const double x = real (z),y = imag (z),l_r = real (l),l_i = imag (l),eta_r = real (eta),eta_i = imag (eta);
		
	if (cwf_real_ptr == 0) cwf_real_ptr = new class Coulomb_wave_functions (is_it_normalized,l_r,eta_r);
	class Coulomb_wave_functions &cwf_real = *cwf_real_ptr;
		
	Complex A_x,dA_x;
	if (is_it_regular) 
		cwf_real.F_dF (x,A_x,dA_x);
	else
		cwf_real.H_dH (1,x,A_x,dA_x);

	const double Bx = real (A_x),dBx = real (dA_x);
	
	double d_l_Bx = 0.0,d_l_dBx = 0.0,d_eta_Bx = 0.0,d_eta_dBx = 0.0;
	if (l_i != 0.0) partial_derivatives (is_it_regular,false,x,d_l_Bx,d_l_dBx);
	if (eta_i != 0.0) partial_derivatives (is_it_regular,true,x,d_eta_Bx,d_eta_dBx);
	
	const double one_over_x = 1.0/x,Coulomb_eq_x = (l_r*(l_r+1)*one_over_x + 2.0*eta_r)*one_over_x - 1.0;
	const double d2Bx = Coulomb_eq_x*Bx;

	B = Complex (Bx,y*dBx + l_i*d_l_Bx + eta_i*d_eta_Bx); 
	dB = Complex (dBx,y*d2Bx + l_i*d_l_dBx + eta_i*d_eta_dBx);
}

void Coulomb_wave_functions::F_dF (Complex const & z,Complex &F,Complex &dF) const
{  
  const double x = real (z),y = imag (z),/*l_r = real (l),*/l_i = imag (l),/*eta_r = real (eta),*/eta_i = imag (eta);

  if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
      && (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
      && (!neg_int_omega_one && !neg_int_omega_minus_one))
    first_order_expansions (true,z,F,dF);
  else if (std::abs (z) <= 0.5)
    F_dF_power_series (z,F,dF);
  else 
  {
    if (real (z) < 0.0) 
      F_dF_with_symmetry_relations (z,F,dF);
    else
    {
      bool is_it_successful = false;
      if (!neg_int_omega_one && !neg_int_omega_minus_one) asymptotic_expansion_F_dF (z,F,dF,is_it_successful);

      if (!is_it_successful) 
      {
	F_dF_direct_integration (z,F,dF,is_it_successful);

	if (!neg_int_omega_one && !neg_int_omega_minus_one && !is_it_successful)
	{
	  const int omega = SIGN(imag (z));
	  const Complex two_I_omega(0.0,2.0*omega),two_I_term = (is_it_normalized) ? (two_I_omega) : (two_I_omega*Cl_eta*Cl_eta);
	  const Complex h_omega = continued_fraction_h (z,omega),h_minus_omega = continued_fraction_h (z,-omega),one_over_two_I_term = 1.0/two_I_term;

	  Complex H_omega,dH_omega;
	  H_dH_direct_integration (omega,z,H_omega,dH_omega);
	  
	  const Complex H_minus_omega = two_I_term/((h_omega - h_minus_omega)*H_omega),dH_minus_omega = h_minus_omega*H_minus_omega;  	  
	  F = (H_omega - H_minus_omega)*one_over_two_I_term;
	  dF = (dH_omega - dH_minus_omega)*one_over_two_I_term;
	}}}
  }

  if (!finite (F) || !finite (dF)) std::cout<<"Numerical failure encountered in F_dF."<<std::endl,exit (1);

  if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0)) F = real (F),dF = real (dF); 
  debut = z,F_debut = F,dF_debut = dF;
}

void Coulomb_wave_functions::G_dG (Complex const & z,Complex &G,Complex &dG) const
{
  const double x = real (z),y = imag (z),/*l_r = real (l),*/l_i = imag (l),/*eta_r = real (eta),*/eta_i = imag (eta);

  if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
      && (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
      && (!neg_int_omega_one && !neg_int_omega_minus_one))
    first_order_expansions (false,z,G,dG);
  else
  {
    if (neg_int_omega_one) 
      H_dH (-1,z,G,dG);
    else if (neg_int_omega_minus_one)
      H_dH (1,z,G,dG);
    else
    {
      Complex F,dF;
      F_dF (z,F,dF);

      Complex H_plus,dH_plus;
      H_dH (1,z,H_plus,dH_plus);

      const Complex I(0.0,1.0);

      if (is_it_normalized)
	G = H_plus - I*F,dG = dH_plus - I*dF;
      else
      {
	const Complex I_Cl_eta_square = I*Cl_eta*Cl_eta;
	
	if ((I_Cl_eta_square == 0.0) || (!finite (I_Cl_eta_square)))
	  G = H_plus - I*exp (2.0*log_Cl_eta + log (F)),dG = dH_plus - I*exp (2.0*log_Cl_eta + log (dF));
	else
	  G = H_plus - I_Cl_eta_square*F,dG = dH_plus - I_Cl_eta_square*dF;
      }
      
      if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0)) G = real (G),dG = real (dG);
    }
  }
}

void Coulomb_wave_functions::H_dH (const int omega,Complex const & z,Complex &H,Complex &dH) const
{
  bool is_it_successful = false;

  Complex H_scaled,dH_scaled;
  if (!neg_int_omega_one && !neg_int_omega_minus_one) asymptotic_expansion_H_dH_scaled (omega,1.0/z,H_scaled,dH_scaled,is_it_successful);

  if (is_it_successful)
  {
    const Complex I_omega(0,omega),log_unscale = I_omega*(z - eta*(M_LN2 + log (z))),unscale = exp (log_unscale);

    if ((unscale == 0.0) || (!finite (unscale)))
      H = exp (log (H_scaled) + log_unscale),dH = exp (log (dH_scaled) + log_unscale);
    else    
      H = H_scaled*unscale,dH = dH_scaled*unscale;
  }
  else
  {
    const double x = real (z),y = imag (z),/*l_r = real (l),*/l_i = imag (l),/*eta_r = real (eta),*/eta_i = imag (eta);
    
    if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
	&& (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
	&& (!neg_int_omega_one && !neg_int_omega_minus_one))
      H_dH_from_first_order_expansions (omega,z,H,dH);
    else
    { 
      //// Replace the following line by : if (!neg_int_omega_one && !neg_int_omega_minus_one) H_dH_with_expansion (omega,z,H,dH,is_it_successful);
      //// if you want H_dH_with_expansion to be used if |l_i| < 1 or |z| > 1.
      if (!neg_int_omega_one && !neg_int_omega_minus_one && (std::abs (l_i) >= 1.0) && (std::abs (z) <= 1.0)) H_dH_with_expansion (omega,z,H,dH,is_it_successful);

      if (!is_it_successful) H_dH_with_F_dF_and_CF (omega,z,H,dH);

      if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0))
      {
	Complex F,dF;
	F_dF (z,F,dF);

	const double omega_norm_functions = (!is_it_normalized) ? (omega*real (Cl_eta)*real (Cl_eta)) : (omega);
	H = Complex (real (H),omega_norm_functions*real (F));
	dH = Complex (real (dH),omega_norm_functions*real (dF));
      }
    }
  }

  if (!finite (H) || !finite (dH)) std::cout<<"Numerical failure encountered in H_dH."<<std::endl,exit (1);
}

void Coulomb_wave_functions::H_dH_scaled (const int omega,Complex const & z,Complex &H_scaled,Complex &dH_scaled) const
{ 
  bool is_it_successful = false;
  if (!neg_int_omega_one && !neg_int_omega_minus_one) asymptotic_expansion_H_dH_scaled (omega,1.0/z,H_scaled,dH_scaled,is_it_successful);

  if (!is_it_successful)
  {
    Complex H,dH;
   
    const double x = real (z),y = imag (z),/*l_r = real (l),*/l_i = imag (l),/*eta_r = real (eta),*/eta_i = imag (eta);

    if (((y != 0.0) || (eta_i != 0.0) || (l_i != 0.0)) 
	&& (std::abs (y) < sqrt_precision*std::min (1.0,x)) && (std::abs (eta_i) < sqrt_precision) && (std::abs (l_i) < sqrt_precision)
	&& (!neg_int_omega_one && !neg_int_omega_minus_one))
      H_dH_from_first_order_expansions (omega,z,H,dH);
    else 
    { 
      //// Replace the following line by : if (!neg_int_omega_one && !neg_int_omega_minus_one) H_dH_with_expansion (omega,z,H,dH,is_it_successful);
      //// if you want H_dH_with_expansion to be used if |l_i| < 1 or |z| > 1.
      if (!neg_int_omega_one && !neg_int_omega_minus_one && (std::abs (l_i) >= 1.0) && (std::abs (z) <= 1.0)) H_dH_with_expansion (omega,z,H,dH,is_it_successful);

      if (!is_it_successful) H_dH_with_F_dF_and_CF (omega,z,H,dH);
      
      if ((y == 0.0) && (eta_i == 0.0) && (l_i == 0.0))
      {  
	Complex F,dF;
	F_dF (z,F,dF);

	const double omega_norm_functions = (!is_it_normalized) ? (omega*real (Cl_eta)*real (Cl_eta)) : (omega);
	H = Complex (real (H),omega_norm_functions*real (F));
	dH = Complex (real (dH),omega_norm_functions*real (dF));
      }
    }
  
    const Complex I_omega(0,omega),log_scale = -I_omega*(z - eta*(M_LN2 + log (z))),scale = exp (log_scale);
    
    if ((scale == 0.0) || (!finite (scale)))
      H_scaled = exp (log (H) + log_scale),dH_scaled = exp (log (dH) + log_scale);
    else
      H_scaled = H*scale,dH_scaled = dH*scale; 
  }

  if (!finite (H_scaled) || !finite (dH_scaled)) std::cout<<"Numerical failure encountered in H_dH_scaled."<<std::endl,exit (1);
}

void Coulomb_wave_functions::F_dF_init (Complex const & z,Complex const & F,Complex const & dF)
{
  debut = z, F_debut = F, dF_debut = dF;
}
