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

#ifndef HEX_CHEBYSHEV
#define HEX_CHEBYSHEV

#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <gsl/gsl_sf.h>

#define sqr(x) gsl_sf_pow_int((x),2)

template <typename Tin, typename Tout> class Chebyshev
{
public:
	
	Chebyshev () {}
	
	/**
	 * Construct Chebyshev 'n'-term approximation of the function 'f'
	 * on the interval (a,b).
	 */
	template <class Functor> Chebyshev (Functor f, int n, Tin a = -1, Tin b = 1)
	{
		N  = n;
		xt = 0.5 * (b + a);
		m  = 0.5 * (b - a);
		
		C.resize(N);
		for (int j = 0; j < N; j++)
		{
			C[j] = 0;
			
			for (int k = 0; k < N; k++)
			{
				double xk = cos(M_PI * (k + 0.5) / N);
				double Tj_xk = cos(M_PI * j * (k + 0.5) / N);
				Tout f_xk = f(unscale(xk));
				
				C[j] += f_xk * Tj_xk;
			}
			
			C[j] *= 2;
			C[j] /= N;
		}
	}
	
	/**
	 * Return full approximation value.
	 */
	Tout operator() (Tin x) const
	{
		Tout ret = 0.5 * C[0];
		double xp = scale (x);
		
		for (int k = 1; k < N; k++)
		{
			double Tk_x = cos(k * acos(xp));
			ret += C[k] * Tk_x;
		}
		return ret;
	}
	
	/**
	 * Use Clenshaw recurrence formula for evaluation of 'm' terms.
	 * The formula has the advantage of not evaluating goniometric funtions.
	 */
	Tout clenshaw (Tin x, int m) const
	{
		Tout d_j = 0, d_jp1 = 0, d_jp2 = 0;
		double one_x = scale(x);
		double two_x = 2 * one_x;
		
		for (int j = m - 1; j >= 1; j--)
		{
			d_j   = two_x * d_jp1 - d_jp2 + C[j];
			
			d_jp2 = d_jp1;
			d_jp1 = d_j;
		}
		
		d_j = one_x * d_jp1 - d_jp2 + 0.5 * C[0];
		
		return d_j;
	}
	
	/**
	 * Get index of the first Chebyshev approximation coefficient
	 * that is smaller than 'eps' times the sum of fabs(C[k])
	 * truncated after this term. If no such term exists, the total
	 * count of terms is returned.
	 * 
	 * NOTE: The Chebyshev polynomial corresponding to the last considered
	 * term can be negligible near some evaluation point 'x'. In that case,
	 * its contribution might be shadowed by the contribution of the following 
	 * polynomial. So, the evaluated result can have worse precision than the
	 * requested 'eps'. Nevertheless, the more polynomials get involved, the
	 * less is this fact problematic.
	 */
	int tail (double eps) const
	{
		double sum = std::abs(0.5 * C[0]);
		double abs_Ck;
		
		for (int k = 1; k < N; k++)
		{
			abs_Ck = std::abs(C[k]);
			sum += abs_Ck;
			
			if (abs_Ck < eps * sum)
				return k + 1;
		}
		
		return N;
	}
	
	/**
	 * Return approximated value where the terms of magnitude less than 'eps'
	 * are discarded.
	 */
	Tout approx (double x, double eps, int * n = 0) const
	{
		Tout ret = 0.5 * C[0];
		double xp = scale(x);
		
		int k;
		for (k = 1; k < N; k++)
		{
			double Tk_x = cos(k * acos(xp));
			Tout delta = C[k] * Tk_x;
			ret += delta;
			
			if (std::abs(delta) / std::abs(ret) < eps)
			{
				k++;
				break;
			}
		}
		
		if (n != 0)
			*n = k;
			
		return ret;
	}
	
	/**
	 * Return Chebyshev aproximation of the function primitive to the
	 * stored Chebyshev approximation.
	 */
	Chebyshev integrate () const
	{
		Chebyshev ret;
		ret.xt = xt;
		ret.m  = m;
		ret.N  = N;
		
		ret.C.resize(ret.N);
		ret.C[0] = 0;
		ret.C[N-1] = ret.m * C[N-2] / (2.*(N-2.));
		
		for (int i = 1; i < N - 1; i++)
			ret.C[i] = ret.m * (C[i-1] - C[i+1]) / (2.*i);
		
		return ret;
	}
	
	/**
	 * Get Chebyshev root in the interval (x1,x2).
	 * \param N Order of the polynomial.
	 * \param k index of the root.
	 * \param x1 Left bound of the interval.
	 * \param x2 Right bound of the interval.
	 */
	static Tin root (int N, int k, Tin x1 = 0., Tin x2 = 1.)
	{
		return x1 + 0.5 * (1. + cos(M_PI * (k + 0.5) / N)) * (x2 - x1);
	}
	
private:
	
	/// map interval (xt-m,xt+m) to (-1,1)
	inline double scale(Tin x) const { return (x - xt) / m; }
	
	/// map interval (-1,1) to (xt-m,xt+m)
	inline Tin unscale(double x) const { return (xt + m*x); }
	
	/// Chebyshev coefficients
	std::vector<Tout> C;
	
	/// coefficient number
	int N;
	
	/// approximation interval center
	Tin xt;
	
	/// approximation interval half-width
	Tin m;
};

/**
 * Limit
 * \f[
 *     \lim_{t \rightarrow x} F(t)
 * \f]
 */
template <class Functor>
auto lim (Functor F, double x, int * n = nullptr) -> decltype(F(0.))
{
	// shorthands for number types
	typedef decltype(F(0.)) FType;
	
	// initial position
	double x0 = (x != 1.) ? 1. : 0.5;
	
	// function value
	FType f0 = F(x0), f;
	
	// main loop
	int i;
	for (i = 0; n == nullptr or i <= *n; i++)
	{
		// advance x0
		if (x == std::numeric_limits<double>::infinity())
			x0 *= 2.;
		else if (x == -std::numeric_limits<double>::infinity())
			x0 *= 2.;
		else
			x0 = 0.5 * (x + x0);
		
		// evaluate function
		f = F(x0);
		
		// check convergence
		if (f == f0)
			break;
		else
			f0 = f;
	}
	
	// return
	if (n != nullptr) *n = i;
	return f;
}

/**
 * Clenshaw-Curtis quadrature for finite interval (a,b).
 * \param F Function to integrate of the signature FType(*)(double x).
 * \param x1 Left bound (allowed infinite).
 * \param x2 Right bound (allowed infinite).
 * \param eps Tolerance.
 * \param n On output, order of the last Chebyshev expansion used.
 */
template <class Functor> 
auto ClenshawCurtis_ff (Functor F, double x1, double x2, double eps = 1e-8, int * n = nullptr) -> decltype(F(0.))
{
	// shorthands for number types
	typedef decltype(F(0.)) FType;
	
	// check interval bounds
	if (x1 == x2)
		return FType(0);
	if (x1 > x2)
		return -ClenshawCurtis_ff(F, x2, x1, eps, n);
	
	// scaled F
	auto f = [&](double x) -> FType {
		return F(x1 + 0.5 * (1.0 + x) * (x2 - x1));
	};
	
	// Chebyshev coefficients
	std::vector<FType> coefs;
	
	// evaluated function f
	std::vector<FType> fvals_prev, fvals;
	
	// integral approximation
	FType sum_prev = std::numeric_limits<double>::quiet_NaN();
	
	// convergence loop
	for (int N = 2; ; N *= 2)
	{
		double pi_over_N = M_PI / N;
		
		fvals.resize(N + 1);
		
		// is this the first iteration?
		if (coefs.empty())
		{
			// evaluate f everywhere
			for (int k = 0; k < N; k++)
			{
				fvals[k] = f(cos(k * pi_over_N));
				
				if (not finite(fvals[k]))
					fvals[k] = 0;
			}
			fvals[N] = f(-1);
			
			if (not finite(fvals[N]))
				fvals[N] = 0;
		}
		else
		{
			// evaluate just the new half, recycle older evaluations
			for (int k = 0; k < N; k++)
			{
				fvals[k] = (k % 2 == 0) ? fvals_prev[k/2] : f(cos(k * pi_over_N));
				
				if (not finite(fvals[k]))
					fvals[k] = 0;
			}
			fvals[N] = fvals_prev[N/2];
		}
		
		coefs.resize(N + 1);
		
		// compute coefficients
		for (int j = 0; j <= N; j++)
		{
			// first and last coef just in half
			if (j % 2 == 0)
				coefs[j] = 0.5 * (fvals[0] + fvals[N]);
			else
				coefs[j] = 0.5 * (fvals[0] - fvals[N]);
			
			// rest of the coefs
			for (int k = 1; k < N; k++)
				coefs[j] += fvals[k] * cos(j * k * pi_over_N);
		}
		
		// sum the quadrature rule
		FType sum = 0.5 * (coefs[0] - coefs[N] / (N*N - 1.));
		for (int twok = 2; twok < N; twok += 2)
			sum -= coefs[twok] / (twok*twok - 1.);
		
		// check convergence
		if (std::abs(sum - FType(2.) * sum_prev) < eps * std::abs(sum))
		{
			if (n != nullptr)
				*n = N;
			return FType(2. * (x2 - x1) / N) * sum;
		}
		
		// save function evaluations and sum
		fvals_prev = fvals;
		sum_prev = sum;
	}
}

/**
 * Clenshaw-Curtis quadrature, main interface.
 * \param F Function to integrate of the signature FType(*)(double x).
 * \param x1 Left bound (allowed infinite).
 * \param x2 Right bound (allowed infinite).
 * \param L Range parameter.
 * \param eps Tolerance.
 * \param n On output, order of the last Chebyshev expansion used.
 */
template <class Functor> 
auto ClenshawCurtis (Functor F, double x1, double x2, double L = 1., double eps = 1e-8, int * n = nullptr) -> decltype(F(0.))
{
	// shorthands for number types
	typedef decltype(F(0.)) FType;
	
	// if both bounds are infinite, call a specialized function
	if (not finite(x1) and not finite(x2))
		return ClenshawCurtis_ii(F, L, eps, n);
	
	// lower bound is infinite
	if (not finite(x1))
	{
		// transform the function
		auto Gwrap = [&](double x) -> FType {
			
			// the functor
			auto G = [&](double t) -> FType {
				double z = L * (1. + t) / (1. - t);
				double w = 2. * L / sqr(1. - t);
				return w * F(x2 - z);
			};
			
			// evaluate
			return (x == 1.) ? lim(G,x) : G(x);
		};
		
		// integrate
		return -ClenshawCurtis_ff(Gwrap, -1., 1., eps, n);	// (-∞,x2)->(1,-1)
	}
	
	// upper bound is infinite
	if (not finite(x2))
	{
		// transform the function
		auto Gwrap = [&](double x) -> FType {
			
			// the functor
			auto G = [&](double t) -> FType {
				double z = L * (1. + t) / (1. - t);
				double w = 2. * L / sqr(1. - t);
				return w * F(x1 + z);
			};
			
			// evaluate
			return (x == 1.) ? lim(G,x) : G(x);
		};
		
		// integrate
		return ClenshawCurtis_ff(Gwrap, -1., 1., eps, n);	// (x1,+∞)->(-1,1)
	}
	
	// both bounds are finite
	return ClenshawCurtis_ff(F, x1, x2, eps, n);
}

/**
 * Clenshaw-Curtis quadrature for infinite-infinite interval (-∞,+∞).
 * \param F Function to integrate of the signature FType(*)(double x).
 * \param L Range parameter.
 * \param eps Tolerance.
 * \param n On output, order of the last Chebyshev expansion used.
 */
template <class Functor> 
auto ClenshawCurtis_ii (Functor F, double L = 1., double eps = 1e-8, int * n = nullptr) -> decltype(F(0.))
{
	// shorthand for function output type
	typedef decltype(F(0.)) FType;
	
	// function values, new and previous
	std::vector<FType> fvals, fvals_prev;
	
	// weights, new and previous
	std::vector<double> weights, weights_prev;
	
	// previous integral
	FType sum_prev = std::numeric_limits<double>::quiet_NaN();
	
	// main loop
	for (int N = 2; ; N *= 2)
	{
		fvals.resize(N);
		weights.resize(N);
		
		// precompute values
		if (fvals_prev.empty())
		{
			// compute all values
			for (int i = 1; i <= N - 1; i++)
			{
				double x = i * M_PI / N;
				fvals[i] = F(L / tan(x));
				if (not finite(fvals[i]))
					fvals[i] = 0;
				weights[i] = 1. / sqr(sin(x));
			}
		}
		else
		{
			// compute new values only
			for (int i = 1; i < N - 1; i++)
			{
				if (i % 2 == 0)
				{
					fvals[i] = fvals_prev[i/2];
					weights[i] = weights_prev[i/2];
				}
				else
				{
					double x = i * M_PI / N;
					fvals[i] = F(L / tan(x));
					if (not finite(fvals[i]))
						fvals[i] = 0;
					weights[i] = 1. / sqr(sin(x));
				}
			}
		}
		
		// evaluate the integral
		FType sum = 0.;
		for (int i = 1; i <= N - 1; i++)
			sum += weights[i] * fvals[i];
		if (std::abs(sum - FType(2.) * sum_prev) < eps * std::abs(sum))
		{
			if (n != nullptr) *n = N;
			return FType(L * M_PI / N) * sum;
		}
		
		// save precomputed values
		sum_prev = sum;
		fvals_prev = fvals;
		weights_prev = weights;
	}
}

#endif
