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
#include <vector>

template <class Functor> class Chebyshev
{
public:
	
	/**
	 * Construct Chebyshev 'n'-term approximation of the function 'f'
	 * on the interval (a,b).
	 */
	Chebyshev (Functor f, int n, double a = -1, double b = 1) : F(f), N(n)
	{
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
				double f_xk = F(unscale(xk));
				
				C[j] += f_xk * Tj_xk;
			}
			
			C[j] *= 2;
			C[j] /= N;
			
			std::cout << "C[" << j << "] = " << C[j] << std::endl;
		}
	}
	
	/**
	 * Return full approximation value.
	 */
	double operator() (double x) const
	{
		double ret = 0.5 * C[0];
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
	double clenshaw (double x, int m) const
	{
		double d_j = 0, d_jp1 = 0, d_jp2 = 0;
		double one_x = scale(x);
		double two_x = 2 * scale(x);
		
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
		double sum = fabs(0.5 * C[0]);
		
		for (int k = 1; k < N; k++)
		{
			double fabs_Ck = fabs(C[k]);
			sum += fabs_Ck;
			
			if (fabs_Ck < eps * sum)
				return k + 1;
		}
		
		return N;
	}
	
	/**
	 * Return approximated value where the terms of magnitude less than 'eps'
	 * are discarded.
	 */
	double approx (double x, double eps, int * n = 0) const
	{
		double ret = 0.5 * C[0];
		double xp = scale(x);
		
		int k;
		for (k = 1; k < N; k++)
		{
			double Tk_x = cos(k * acos(xp));
			double delta = C[k] * Tk_x;
			ret += delta;
			
			if (fabs(delta) / fabs(ret) < eps)
			{
				k++;
				break;
			}
		}
		
		if (n != 0)
			*n = k;
			
		return ret;
	}
	
private:
	
	/// initialize
	void init(double a, double b);
	
	/// map interval (xt-m,xt+m) to (-1,1)
	inline double scale(double x) const { return (x - xt) / m; }
	
	/// map interval (-1,1) to (xt-m,xt+m)
	inline double unscale(double x) const { return (xt + m*x); }
	
	/// approximated function
	Functor F;
	
	/// Chebyshev coefficients
	std::vector<double> C;
	
	/// coefficient number
	int N;
	
	/// approximation interval center
	double xt;
	
	/// approximation interval half-width
	double m;
};

#endif
