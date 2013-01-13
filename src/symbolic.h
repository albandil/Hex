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

#ifndef HEX_SYMBOLIC
#define HEX_SYMBOLIC

#include <vector>
#include <cln/cln.h>

#define GF_NONE	0
#define GF_SIN	1
#define GF_COS	2

/**
 * \brief Symbolic term structure
 * 
 * Contains information about one term of the form
 * \f[
 *       k x^a \mathrm{gon}\,(bx) \exp(-cx)
 * \f]
 */
class SymbolicTerm
{
public:
	
	//
	// constructors
	//
	
	SymbolicTerm() : ki(1), kr(0), a(0), gf(GF_NONE), b(0), c(0) {}
	SymbolicTerm(cln::cl_RA x) : ki(1), kr(x), a(0), gf(GF_NONE), b(0), c(0) {}
	
	//
	// data fields
	//
	
	/// constant irrational multiplication factor
	double ki;
	
	/// constant rational multiplication factor
	cln::cl_RA kr;
	
	/// exponent
	int a;
	
	/// goniometric function (none, sin or cos)
	int gf;
	
	/// goniometric function wave number
	double b;
	
	/// exponential factor
	cln::cl_RA c;
	
};

SymbolicTerm operator + (SymbolicTerm const & A, SymbolicTerm const & B);

/**
 * \brief Sum of symbolic terms.
 */
class SymbolicPoly
{
public:
	
	//
	// constructors
	//
	
	SymbolicPoly() {}
	SymbolicPoly(int n) : terms(n) {}
	SymbolicPoly(SymbolicTerm const & T) : terms({T}) {}
	SymbolicPoly(SymbolicPoly const & P) : terms(P.terms) {}
	
	//
	// data storage
	//
	std::vector<SymbolicTerm> terms;
	void optimize();
	
	//
	// std::vector-line interface
	//
	
	typedef std::vector<SymbolicTerm>::iterator       iterator;
	typedef std::vector<SymbolicTerm>::const_iterator const_iterator;
	inline iterator begin() { return terms.begin(); }
	inline iterator end()   { return terms.end();   }
	inline const_iterator begin() const { return terms.begin(); }
	inline const_iterator end()   const { return terms.end();   }
	inline SymbolicTerm & front() { return terms.front(); }
	inline SymbolicTerm & back()  { return terms.back();  }
	inline SymbolicTerm const & front() const { return terms.front(); }
	inline SymbolicTerm const & back()  const { return terms.back();  }
	inline void push_back(SymbolicTerm const & T) { terms.push_back(T); }
	inline void insert(iterator a, const_iterator u, const_iterator v) { terms.insert(a,u,v); }
	inline size_t size() const { return terms.size(); }
	inline SymbolicTerm & operator [] (size_t i) { return terms[i]; }
	inline SymbolicTerm const & operator [] (size_t i) const { return terms[i]; }
};

/// arithmetic
SymbolicPoly operator + (SymbolicPoly const & P, SymbolicPoly const & Q);
SymbolicPoly operator - (SymbolicPoly const & P, SymbolicPoly const & Q);
SymbolicPoly operator * (SymbolicPoly const & P, SymbolicPoly const & Q);

/// Laguerre polynomial
SymbolicPoly Laguerre(int k, int s);

/// Hydrogen radial function normalization factor
SymbolicTerm HydrogenN(int n, int l);

/// Hydrogen radial function.
SymbolicPoly HydrogenP(int n, int l);

/// Hydrogen sturmian function.
SymbolicPoly HydrogenS(int n, int l, cln::cl_RA lambda);

/// Riccati-Bessel function
SymbolicPoly RiccatiBessel(int l, double k);

/// integrals
SymbolicPoly integrate_inf(SymbolicPoly const & P);
SymbolicPoly integrate_low(SymbolicPoly const & P);
SymbolicTerm integrate_full(SymbolicPoly const & P);

/// evaluate the polynomial
double eval(SymbolicPoly const & P, double r);

/// formatted output
void write(SymbolicPoly const & p);

#endif
