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

#ifndef HEX_FORBIDDEN_WAVE
#define HEX_FORBIDDEN_WAVE

#include "complex.h"
#include "potential.h"
#include "ode.h"
#include "specf.h"

/**
 * \brief Irregular forbidden wave information.
 * 
 */
class ForbiddenWave : public RadialFunction<double>
{
public:
    
    // constructors
    // @{
    ForbiddenWave(double _kn, int _ln, DistortingPotential const & _U);
    ForbiddenWave(ForbiddenWave const& W) { *this = W; }
    // @}
    
    ~ForbiddenWave();
    
    // assignment
    ForbiddenWave operator= (ForbiddenWave const& W);
    
    /**
     * Evaluate distorted wave.
     */
    double operator()(double x) const;
    
    /// Classical turning point.
    double getTurningPoint () const { return r0; }
    
    /// Near-zero asymptotic behaviour.
    std::pair<double,int> getZeroAsymptotic (double x) const;
    
    /**
     * \brief Return derivatives from the distorted wave equation.
     * 
     */
    int derivs(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
    int derivs0(double x, size_t nv, const o2scl::ovector_base& y, o2scl::ovector_base& dydx);
    
    /**
     * Export data to file using \ref write_array.
     */
    void toFile(const char * filename) const;
    
    void scale(bool s);
    
    /// (debuging parameter) number of evaluations
    mutable unsigned Evaluations;
    
private:
    
    /// distorting potential
    DistortingPotential U;
    
    /// interpolator
    gsl_interp *interpolator, *interpolator0;
    
    /// wavenumber of the distorted wave
    double kn;
    
    /// angular momentum of the distorted wave
    int ln;
    
    /// sample count
    int samples, samples0;
    
    /// discretization step
    double h, h0;
    
    /// grid
    rArray grid, grid0;
    
    /// samples
    rArray array, array0;
    
    bool Scaled;
    
    /// classical turning point, far radius
    double r0, rf;
};

#endif
