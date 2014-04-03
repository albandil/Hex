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

#ifndef HEX_PWBA2_MULTIPOT
#define HEX_PWBA2_MULTIPOT

/**
 * @brief Multipole potential.
 * 
 * Defined as
 * @f[
 *     V_{\alpha \beta}^\lambda (r) = \int_0^\infty
 *     \left<\psi_\alpha(r')\right|
 *     \left(\frac{r_<^\lambda}{r_>^{\lambda+1}} - \frac{1}{r}\right)
 *     \left|\psi_\beta(r')\right> \,
 *     \mathrm{d}r' \ .
 * @f]
 * For bound-bound potential the result is
 * @f[
 *     \sum_{k=1}^{N_\alpha+L_\alpha} \sum_{l=1}^{N_\beta+L_\beta}
 *     a_k b_l \sum_{m=0}^{L_\alpha+L_\alpha+k+l}
 *     \frac{1}{m!}\left(\frac{1}{N_\alpha}+\frac{1}{N_\beta}\right)^m \ .
 * @f]
 */
class MultipolePotential
{
    public:
        
        // constructors
        MultipolePotential (int lambda, int Na,    int La, int Nb,    int Lb);
        MultipolePotential (int lambda, double Ka, int La, int Nb,    int Lb);
        MultipolePotential (int lambda, int Na,    int La, double Kb, int Lb);
        
        // evaluate the potential
        double operator() (double x);
        
        // different types of the potential
        typedef enum
        {
            bound_bound,
            free_bound,
            bound_free
        }
        type;
        
    private:
        
        // multipole moment (transferred angular momentum)
        int lambda;
        
        // final or initial bound state quantum numbers
        int Na, La, Nb, Lb;
        
        // final or initial free state wavenumbers
        double Ka, Kb;
};

#endif
