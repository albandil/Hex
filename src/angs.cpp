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
#include <map>
#include <tuple>
#include <gsl/gsl_sf.h>

#define Wigner3j(a,b,c,d,e,f) gsl_sf_coupling_3j(2*(a),2*(b),2*(c),2*(d),2*(e),2*(f))
#define Wigner6j(a,b,c,d,e,f) gsl_sf_coupling_6j(2*(a),2*(b),2*(c),2*(d),2*(e),2*(f))

double computef(int lambda, int l1, int l2, int l1p, int l2p, int L)
{
	return pow(-1, L + l2 + l2p) * sqrt((2*l1 + 1) * (2*l2 + 1) * (2*l1p + 1) * (2*l2p + 1))
				* Wigner6j(l1, l2, L, l2p, l1p, lambda)
				* Wigner3j(l1, lambda, l1p, 0, 0, 0)
				* Wigner3j(l2, lambda, l2p, 0, 0, 0);
}

std::map<long double, long double> wdfact;

// inline long double dfact(long double __x)
// {
// 	if (wdfact.find(__x) == wdfact.end())
// 	{
// 		if((__x < 0.00001) && (__x >= 0.0)) return 1.;
// 		if(__x < 0) return 0.;
// 		return wdfact[__x] = __x*dfact(__x - 1.);
// 	}
// 	else
// 	{
// 		return wdfact[__x];
// 	}
// }

inline long double dfact(long double x)
{
	if(x < 0)
		return 0.;
	
	if(x < 0.00001)	// x == 0 (+ rounding errors)
		return 1.;
	
	return x * dfact(x - 1.);
}

double ClebschGordan(int __j1, int __m1, int __j2, int __m2, int __J, int __M)
{
	if((__m1 + __m2) != __M) return 0.;
	if(abs(__m1) > __j1) return 0;
	if(abs(__m2) > __j2) return 0;

	// convert to pure integers (each 2*spin)
	int j1 = (int)(2.*__j1);
	int m1 = (int)(2.*__m1);
	int j2 = (int)(2.*__j2);
	int m2 = (int)(2.*__m2);
	int J = (int)(2.*__J);
	int M = (int)(2.*__M);

	long double n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,A,exp;
	int nu = 0;
  
	long double sum = 0;
	while (((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 ))
		nu++;
	while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0) && ((n1=(j2+J+m1)/2-nu) >= 0 ))
	{
		d3=((j1-j2-M)/2+nu);
		n2=((j1-m1)/2+nu);
		d0=dfact((double) nu);
		exp=nu+(j2+m2)/2;
		n0 = (double) pow(-1.,exp);
		sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
		nu++;
	}

	if (sum == 0) return 0;

	n0 = J+1;
	n1 = dfact((double) (J+j1-j2)/2);
	n2 = dfact((double) (J-j1+j2)/2);
	n3 = dfact((double) (j1+j2-J)/2);
	n4 = dfact((double) (J+M)/2);
	n5 = dfact((J-M)/2);
  
	d0 = dfact((double) (j1+j2+J)/2+1);
	d1 = dfact((double) (j1-m1)/2);
	d2 = dfact((double) (j1+m1)/2);
	d3 = dfact((double) (j2-m2)/2);
	d4 = dfact((double) (j2+m2)/2);
  
	A = ((long double) (n0*n1*n2*n3*n4*n5))/((long double) (d0*d1*d2*d3*d4));
	
	return sqrtl(A)*sum;           
}

double Gaunt(int l1, int m1, int l2, int m2, int l, int m)
{
	// dictionary
	static std::map<std::tuple<int,int,int,int,int,int>,double> dict;
	
	// dictionary key
	std::tuple<int,int,int,int,int,int> key = std::make_tuple(l1,m1,l2,m2,l,m);
	
	// try to find this Gaunt's coefficient in the dictionary
	if (dict.find(key) != dict.end())
		return dict[key];
	
	// compute the value and store it to the dictionary
	dict[key] = sqrt((2*l1+1)*(2*l2+1) / (4*M_PI*(2*l+1))) *
		ClebschGordan(l1, m1, l2, m2, l, m) *
		ClebschGordan(l1,  0, l2,  0, l, 0);
	return dict[key];
}

int triangle_count(int L, int maxl)
{
	int n = 0;
	
	for (int l1 = 0; l1 <= maxl; l1++)
		for (int l2 = 0; l2 <= maxl; l2++)
			if (std::abs(l1-l2) <= L and l1+l2 >= L)
				n++;
	
	return n;
}
