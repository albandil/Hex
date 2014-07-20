/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_DIOPHANTINE
#define HEX_DIOPHANTINE

#include <cmath>
#include <map>

#include "arrays.h"

namespace dioph
{

/**
 * @brief Diophantine quadrature modes.
 * 
 * List of available quadrature modes. The name pattern is
 * <pre>
 * d&lt;dimension&gt;n&lt;order&gt;
 * </pre>
 * The 'dimension' used has to be compatible with the function that is to be
 * integrated. The 'order' is the number of evaluation points.
 */
typedef enum
{
    // dim = 2
    d2n52, d2n538, d2n1154, d2n3722, d2n6044,
    // dim = 3
    d3n52, d3n538, d3n1154, d3n3722, d3n6044,
    // dim = 4
    d4n1154, d4n3722, d4n6044,
    // dim = 5
    d5n1154, d5n3722, d5n6044,
    // dim = 6
    d6n1154, d6n2008, d6n3722, d6n6044, d6n9644,
    // dim = 7
    d7n1154, d7n3722, d7n6044,
    // dim = 8
    d8n3708, d8n6044, d8n9644,
    // dim = 9
    d9n3722, d9n6044, d9n9644,
    // dim = 10
    d10n3722, d10n6044, d10n9644,
    // dim = 11
    d11n3722, d11n6044, d11n9644,
    // dim = 12
    d12n3722, d12n6044, d12n9644,
}
id;

/// Internal list of evaluation nodes.
extern const std::map<id,iArray> nodes;

} // end of namespace "dioph"

/**
 * @brief Diophantine quadrature.
 * 
 * This routine will integrate an n-dimensional function of a real variable
 * (but the return type can be arbitrary number type) using the Diophantine
 * method presented in
 * <center><pre>
 * Conroy, H., J. Chem. Phys. 47 (12) (1967), 5307-5318.
 * </pre></center>
 * The integration domain is the unit cube (0,1)^n in ℝ^n. It is assumed that
 * the supplied function is symmetrical with respect to reflection,
 * @f[
 *    F(x_1, \dots, x_i, \dots, x_n) = F(n_1, \dots, -x_i, \dots, x_n) \ ,
 * @f]
 * which can be easily secured by using @f$ F @f$ that depends only on the
 * absolute values of its arguments. In the program implementation the function
 * F should accept two parameters: the dimension and pointer to the list of
 * doubles of length equal to the dimension (the coordinates).
 * 
 * @param F The n-dimensional function to integrate, compatible with the signature
 *           T (*) (int, double const *).
 * @param id Mode of integration (contains dimension and order of quadrature).
 */
template <class T, class Function> T diophantine (Function F, dioph::id id)
{
    // check that this "id" is implemented
    if (dioph::nodes.find(id) == dioph::nodes.end())
        throw exception ("The requested Diophantine id = %d is not implemented.", id);
    
    // get dimension of the quadrature
    int dim = dioph::nodes.at(id)[0];
    
    // get number of samples
    int M = dioph::nodes.at(id)[1];
    
    // get sample points
    double a[dim];
    for (int d = 0; d < dim; d++)
        a[d] = dioph::nodes.at(id)[2+d] / double(M);
    
    // evaluate the function at sample points
    T suma = 0;
    double intpart = 0, x[dim];
    for (int two_m = 1; two_m <= M; two_m += 2)
    {
        // compute evaluation vector
        for (int d = 0; d < dim; d++)
            x[d] = std::abs(2 * std::modf(0.5 * (two_m * a[d] + 1), &intpart) - 1);
        
        // evaluate function at 'x'
        suma += F(dim,x);
    }
    
    // normalize and return the sum
    return (T(2) * suma) / T(M);
}

#endif
