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

#include <fftw3.h>
#include <gsl/gsl_sf.h>

#include "compact.h"

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

        // evaluate nodes and function
        std::vector<double> x(N);
        std::vector<Tout> fval(N);
        for (int k = 0; k < N; k++)
        {
            x[k] = cos(M_PI * (k + 0.5) / N);
            fval[k] = f(unscale(x[k]));
        }

        // compute the coefficients
        C.resize(N);
        for (int j = 0; j < N; j++)
        {
            C[j] = 0;

            for (int k = 0; k < N; k++)
            {
                double Tj_xk = cos(M_PI * j * (k + 0.5) / N);
                C[j] += fval[k] * Tj_xk;
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

            if (abs_Ck <= eps * sum)
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

            if (std::abs(delta) <= eps * std::abs(ret))
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

    /**
     * Write out the coefficients.
     */
    std::string str() const
    {
        std::ostringstream out;

        out << "[";
        if (N > 0)
        {
            for (int i = 0; i < N - 1; i++)
                out << C[i] << ", ";
            out << C[N-1];
        }
        out << "]" << std::endl;

        return out.str();
    }

private:

    /// map interval (xt-m,xt+m) to (-1,1)
    inline double scale(Tin x) const {
        return (x - xt) / m;
    }

    /// map interval (-1,1) to (xt-m,xt+m)
    inline Tin unscale(double x) const {
        return (xt + m*x);
    }

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
 * Clenshaw-Curtis quadrature for finite interval (a,b).
 * \param F Function to integrate of the signature FType(*)(double x).
 * \param x1 Left bound (allowed infinite).
 * \param x2 Right bound (allowed infinite).
 * \param eps Tolerance.
 * \param n On input, maximal subdivision for a single bisection. Use *n = -1
 *          to disable recurrence. On output, evaluations needed for a converged result.
 */
template <class Functor, typename FType>
FType ClenshawCurtis_ff (Functor const & F, double x1, double x2, double eps = 1e-8, int * n = nullptr)
{
    // check interval bounds
    if (x1 == x2)
        return FType(0);
    if (x1 > x2)
        return -ClenshawCurtis_ff<decltype(F),FType>(F, x2, x1, eps, n);

// 	std::cout << "[ClenshawCurtis_ff] x1 = " << x1 << ", x2 = " << x2 << "\n";
	
    // scaled F
    auto f = [&](double x) -> FType {
        return F(x1 + 0.5 * (1.0 + x) * (x2 - x1));
    };

    // Chebyshev coefficients
    std::vector<FType> coefs;

    // evaluated function f
    std::vector<Complex> fvals_prev, fvals;

    // integral approximation
    FType sum_prev = std::numeric_limits<double>::quiet_NaN();

	// subdivision limit
	int MaxN = (n == nullptr) ? 32 : *n;
	
    // convergence loop
    for (int N = 4; MaxN < 0 or N <= MaxN; N *= 2)
    {
        double pi_over_N = M_PI / N;

        fvals.resize(2*N + 1);

        // is this the first iteration?
        if (coefs.empty())
        {
            // evaluate f everywhere
            for (int k = 0; k < N; k++)
            {
                fvals[k] = fvals[2*N-k] = f(cos(k * pi_over_N));

                if (not finite(std::abs(fvals[k])))
                    fvals[k] = fvals[2*N-k] = 0;
            }
            fvals[N] = f(-1);

            if (not finite(std::abs(fvals[N])))
                fvals[N] = 0;
        }
        else
        {
            // evaluate just the new half, recycle older evaluations
            for (int k = 0; k < N; k++)
            {
                fvals[k] = fvals[2*N-k] = (k % 2 == 0) ? fvals_prev[k/2] : f(cos(k * pi_over_N));
				
                if (not finite(std::abs(fvals[k])))
                    fvals[k] = fvals[2*N-k] = 0;
            }
            fvals[N] = fvals_prev[N/2];
        }

        coefs.resize(N + 1);

//         std::cout << "[ClenshawCurtis_ff] eval " << N <<  " coefs..." << std::endl;

#if 0
		// compute coefficients
        for (int j = 0; j <= N; j++)
        {
            // first and last coef just in half
            if (j % 2 == 0)
                coefs[j] = 0.5 * (fvals[0] + fvals[N]).real();
            else
                coefs[j] = 0.5 * (fvals[0] - fvals[N]).real();

            // rest of the coefs
            for (int k = 1; k < N; k++)
                coefs[j] += fvals[k].real() * cos(j * k * pi_over_N);
        }
#else
		// compute coefficients using FFT
		std::vector<Complex> ftraf(2*N);
		fftw_plan plan = fftw_plan_dft_1d (
			2*N,
			reinterpret_cast<fftw_complex*>(&fvals[0]),
			reinterpret_cast<fftw_complex*>(&ftraf[0]),
            FFTW_BACKWARD,
			0
		);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		
		// create type-correct pointer
		FType const * ftraf_ptr = reinterpret_cast<FType*>(&ftraf[0]);
		
		// copy result
		for (int i = 0; i <= N; i++)
		{
			if (typeid(FType) == typeid(Complex))
			{
				// copy whole complex number
 				coefs[i] = 0.5 * (*(ftraf_ptr + i));
			}
			else if (typeid(FType) == typeid(double))
			{
				// copy just real part
 				coefs[i] = 0.5 * (*(ftraf_ptr + 2*i));
			}
			else
			{
				std::cerr << "[ClenshawCurtis_ff] Can't handle datatype.\n";
				throw;
			}
		}
#endif
		

//         std::cout << "[ClenshawCurtis_ff] coefs: [";
//         for (int i = 0; i < N; i++)
//             std::cout << coefs[i] << ",";
//         std::cout << coefs[N] << "]\n";

//         std::cout << "[ClenshawCurtis_ff] sum...\n";

        // sum the quadrature rule
        FType sum = 0.5 * (coefs[0] - coefs[N] / (N*N - 1.));
        for (int twok = 2; twok < N; twok += 2)
            sum -= coefs[twok] / (twok*twok - 1.);

//         std::cout << "[ClenshawCurtis_ff] check convergence...\n";

        // check convergence
        if (std::abs(sum - FType(2.) * sum_prev) <= eps * std::abs(sum))
        {
// 			std::cout << "[ClenshawCurtis_ff] OK, sum = " << sum << std::endl;
            if (n != nullptr)
                *n = N;
			
// 			std::cout << "[ClenshawCurtis_ff] " << N << " " << sum << " " << 2.*sum_prev << " -> " << FType(2. * (x2 - x1) / N) * sum << " \n";
			
            return FType(2. * (x2 - x1) / N) * sum;
        }
        else if (finite(std::abs(sum)) and finite(std::abs(sum_prev)) and std::max(std::abs(sum), std::abs(sum_prev)) <= 1e-15)
		{
			// numerical niose only, return clean zero
// 			std::cout << "[ClenshawCurtis_ff] Return zero.\n";

// 			std::cout << "[ClenshawCurtis_ff] " << N << " " << sum << " " << 2.*sum_prev << " -> 0 \n";

			return FType(0.);
		}
		else
        {
// 			std::cout << "[ClenshawCurtis_ff] " << N << " " << sum << " " << 2.*sum_prev << " -> recycle \n";
        }

        // save function evaluations and sum
        fvals_prev = fvals;
        sum_prev = sum;
    }
    
    // no convergence? -> bisect
//     std::cout << "[ClenshawCurtis_ff] x1 = " << x1 << ", x2 = " << x2 << " -> bisect \n";
    int n1, n2;
	n1 = n2 = (n == nullptr) ? 16 : *n;
    FType i1 = ClenshawCurtis_ff<Functor,FType>(F, x1, (x2+x1)/2, eps, &n1);
	FType i2 = ClenshawCurtis_ff<Functor,FType>(F, (x2+x1)/2, x2, eps, &n2);
	if (n != nullptr) *n = n1 + n2;
// 	std::cout << "[ClenshawCurtis_ff] x1 = " << x1 << ", x2 = " << x2 << " -> recurrent result: " << i1 + i2 << "\n";
	return i1 + i2;
}

/**
 * Clenshaw-Curtis quadrature for infinite-infinite interval (-∞,+∞).
 * \param F Function to integrate of the signature FType(*)(double x).
 * \param L Range parameter.
 * \param eps Tolerance.
 * \param n On output, evaluations needed for 
 *          converged result.
 */
template <class Functor, typename FType>
FType ClenshawCurtis_ii (Functor const & F, double L = 1., double eps = 1e-8, int * n = nullptr)
{
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
                if (not finite(std::abs(fvals[i])))
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
                    if (not finite(std::abs(fvals[i])))
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

/**
 * Clenshaw-Curtis quadrature, main interface.
 * \param F Function to integrate of the signature FType(*)(double x).
 * \param x1 Left bound (allowed infinite).
 * \param x2 Right bound (allowed infinite).
 * \param limit Whether to use limit (limit = true) when evaluating F(x1) and
 *              F(x2) for improper arguments x1=-∞ and/or x2=+∞. Otherwise 
 *              (limit = false) the functor will be simply evaluated at
 *              possibly infinite boundaries and it has to cope itself with the
 *              input.
 * \param L Range parameter.
 * \param eps Tolerance.
 * \param n On input, maximal subdivision for a single bisection. (No effect
 *          for double infinite intervals.) On output, evaluations needed for 
 *          converged result.
 */
template <class Functor, typename FType>
FType ClenshawCurtis (
	Functor const & F, double x1, double x2, 
	bool limit = true, double L = 1., double eps = 1e-8, int * n = nullptr
){
	if (x1 == x2)
		return 0.;
	
    // if both bounds are infinite, call a specialized function
    if (not finite(x1) and not finite(x2))
        return ClenshawCurtis_ii<decltype(F),FType>(F, L, eps, n);

    // lower bound is infinite
    if (not finite(x1))
    {
        // the compactified functor
        CompactIntegrand<decltype(F),FType> G(F, x1, x2, limit, L);

        // integrate
        return -ClenshawCurtis_ff<decltype(G),FType>(G, -1., 1., eps, n);	// (-∞,x2)->(1,-1)
    }

    // upper bound is infinite
    if (not finite(x2))
    {
        // the compactified functor
        CompactIntegrand<decltype(F),FType> G(F, x1, x2, limit, L);

        // integrate
        FType f = ClenshawCurtis_ff<decltype(G),FType>(G, -1., 1., eps, n);	// (x1,+∞)->(-1,1)
//         std::cout << "[ClenshawCurtis] f = " << f << ", returning.\n";
        return f;
    }

    // both bounds are finite
    return ClenshawCurtis_ff<decltype(F),FType>(F, x1, x2, eps, n);
}

#endif
