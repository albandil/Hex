//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_DISTORTED_WAVE
#define HEX_DISTORTED_WAVE

#include <gsl/gsl_interp.h>

#include "complex.h"
#include "potential.h"
#include "special.h"

/**
 * \brief Distorted wave information.
 * 
 * The distorted wave \f$ \phi_{l_n}(k_n, r) \f$ is a solution of the equation
 * \f[
 *     \left(-\frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l_n(l_n+1)}{r^2} + 2U_n(r) \right)
 *     \phi_{l_n} (k_n, r) = k_n^2 \phi_{l_n} (k_n, r) 
 * \f]
 * which is regular at origin and which satisfies the boundary condition
 * \f[
 *     \phi_{l_n} (k_n, r) \propto \hat{j}_{l_n}(k_n r) + t_{l_n}(k_n) \mathrm{i}\hat{h}_{l_n}^{(+)}(k_n r) \ .
 * \f]
 * Most of the work is being done in the constructor. The radial interval between 
 * the origin and the far radius of DistortingPotential is divided into an equidistant
 * grid (with spacing adjusted to the energy of projectile),
 * on which the solution is sought. Integration starts from \f$ r = h \f$, where \f$ h \f$
 * is a small number compared to grid size. The initial conditions here are chosen as
 * \f[
 *     \phi_{l_n}(k_n, h) \simeq h \cdot 10^{-5} \ ,
 * \f]
 * \f[
 *     \phi_{l_n}'(k_n, h) \simeq (l_n + 1) \cdot 10^{-5} \ ,
 * \f]
 * which is in accord with the behaviour of Ricatti-Bessel functions near origin. Resulting
 * solution is real and using the standard formula
 * \f[
 *     \tan \delta_{l_n} = \frac{k_n \cos(k_n R - l_n\pi/2) - D \sin(k_n R - l_n\pi/2)}
 *                                {k_n \sin(k_n R - l_n\pi/2) + D \cos(k_n R - l_n\pi/2)}
 * \f]
 * the phase shift is determined. This number sets the phase factor \f$ \mathrm{e}^{\mathrm{i}\delta_{l_n}} \f$
 * of the whole wave.
 */
class DistortedWave : public special::RadialFunction<double>
{
public:
    
    // constructors
    // @{
    DistortedWave(double _kn, int _ln, DistortingPotential const & _U);
    DistortedWave(DistortedWave const& W) { *this = W; }
    // @}
    
    // destructor
    ~DistortedWave();
    
    // assignment
    DistortedWave operator= (DistortedWave const& W);
    
    /**
     * Evaluate distorted wave.
     */
    double operator() (double x) const;
    
    /**
     * Return the phase factor \f$ \mathrm{e}^{\mathrm{i}\delta_{l_n}} \f$.
     */
    Complex getPhasef() const { return Complex(cos(phase),sin(phase)); }
    
    /**
     * Return the phase \f$ \delta_{l_n} \f$.
     */
    double getPhase() const { return phase; }
    
    /// Wavenumber.
    double k() const { return kn; }
    
    /// Angular momentum.
    int l() const { return ln; }
    
    /// Classical turning point.
    double getTurningPoint () const { return r0; }
    
    /// Near-zero asymptotic behaviour.
    std::pair<double,int> getZeroAsymptotic (double x) const;
    
    /// Export data to file using \ref write_array.
    void toFile(const char * filename) const { write_array(grid, array, filename); }
    
    /**
     * Return radius from which the asymptotic form is used.
     */
    double farRadius() const { return rf; }
    
    std::size_t sampleCount() const { return grid.size(); }
    
    /// (debuging parameter) number of evaluations
    mutable unsigned Evaluations;
    
private:
    
    // distorting potential
    DistortingPotential U;
    
    // interpolator
    gsl_interp *interpolator, *interpolator0;
    
    // distorted wave input parameters
    double kn;      // wavenumber of the distorted wave
    int ln;         // angular momentum of the distorted wave
    
    // classical turning point and far radius
    double r0, rf;
    
    // distorted wave computed attributes
    int samples, samples0;  // sample count
    double h, h0;           // discretization step
    
    // distorted wave data
    double phase;           // phase shift
    rArray grid, grid0;     // grid
    rArray array, array0;   // samples
};

#endif
