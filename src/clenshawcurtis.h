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

#ifndef HEX_CLENSHAWCURTIS
#define HEX_CLENSHAWCURTIS

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include <fftw3.h>
#include <gsl/gsl_sf.h>

#include "arrays.h"
#include "compact.h"
#include "misc.h"

/**
 * @brief Clenshaw-Curtis quadrature.
 * 
 * Clenshaw-Curtis integrator class. Constructor of the class accepts the
 * function that will be integrated. The actual quadrature is done on the
 * call to ClenshawCurtis::integrate. Clenshaw-Curtis quadrature evaluates
 * the function in the Chebyshev nodes and uses sophisticated algorithm to
 * extract the value of the quadrature. For more details see Numerical
 * recipes (3rd ed.). The present implementation uses fast Fourier transform
 * from FFTW for fast computation.
 */
template <class Functor, typename FType> class ClenshawCurtis
{
public:
    
    /**
     * Constructor.
     * @param  f Function to integrate of the signature FType(*)(double x).
     */ 
    ClenshawCurtis (Functor const & f) : F(f), EpsRel(1e-8), EpsAbs(1e-12), 
        Limit(false), Recurrence(true), NNest(5), NStack(5), L(1.0), Verbose(false),
        vName("[ClenshawCurtis_ff]"), Log(std::cout.rdbuf()), Throw(true) {}
    
    /// Get relative tolerance.
    inline double eps() const { return EpsRel; }
    
    /// Set relative tolerance.
    inline void setEps(double epsrel) { EpsRel = epsrel; }
    
    /// Get absolute tolerance.
    inline double tol() const { return EpsAbs; }
    
    /// Set absolute tolerance.
    inline void setTol(double epsabs) { EpsAbs = epsabs; }
    
    /// Get compacification parameter.
    inline double getRange() const { return L; }
    
    /// Set compacification parameter.
    inline void setRange(double l) { L = l; }
    
    /// Get 2-log of maximal subdivision interval count.
    inline int subdiv() const { return NNest; }
    
    /// Set 2-log of maximal subdivision interval count.
    inline void setSubdiv(int nlevel) { NNest = nlevel; }
    
    /// Get subdivision level limit.
    inline int stack() const { return NStack; }
    
    /// Set subdivision level limit.
    inline void setStack(int nlevel) { NStack = nlevel; }
    
    /// Get limit flag.
    inline bool lim() const { return Limit; }
    
    /// Set limit flag.
    inline void setLim(bool limit) { Limit = limit; }
    
    /// Get recurrence flag.
    inline bool Rec() const { return Recurrence; }
    
    /// Set recurrence flag.
    inline void setRec(bool recurrence) { Recurrence = recurrence; }
    
    /// Get verbose flag.
    inline bool verbose() const { return Verbose; }
    
    /// Set verbose flag.
    inline void setVerbose
    (
        bool verbose,
        std::string name = "ClenshawCurtis_ff",
        std::ostream & stream = std::cout
    )
    {
        Verbose = verbose;
        vName = name;
        Log.rdbuf(stream.rdbuf());
    }
    
    /// Set warn flag.
    inline void setThrowAll(bool t) { Throw = t; }
    
    /// Get warn flag.
    inline bool throwall() const { return Throw; }
    
    /**
     * @brief General integration.
     * 
     * Clenshaw-Curtis quadrature, main interface.
     * 
     * @param x1 Left bound (allowed infinite).
     * @param x2 Right bound (allowed infinite).
     * @param n On input, maximal subdivision for a single bisection. (No effect
     *          for double infinite intervals.) On output, evaluations needed for 
     *          converged result.
     */
    FType integrate (double x1, double x2, int * n = nullptr) const
    {
        if (x1 == x2)
            return 0.;
        
        // if both bounds are infinite, call a specialized function
        if (not std::isfinite(x1) and not std::isfinite(x2))
            return integrate_ii(n);

        // lower bound is infinite
        if (not std::isfinite(x1))
        {
            // the compactified functor
            CompactIntegrand<decltype(F),FType> G(F, x1, x2, Limit, L);
            
            // quadrature system
            ClenshawCurtis<decltype(G),FType> cc_G(G);
            cc_G.setEps(EpsRel);
            cc_G.setTol(EpsAbs);
            cc_G.setLim(Limit);
            cc_G.setRec(Recurrence);
            cc_G.setSubdiv(NNest);
            cc_G.setStack(NStack);
            cc_G.setRange(L);
            cc_G.setVerbose(Verbose, vName, Log);
            
            // integrate
            return -cc_G.integrate_ff(-1., 1., n);    // (-∞,x2)->(1,-1)
        }

        // upper bound is infinite
        if (not std::isfinite(x2))
        {
            // the compactified functor
            CompactIntegrand<decltype(F),FType> G(F, x1, x2, Limit, L);
            
            // quadrature system
            ClenshawCurtis<decltype(G),FType> cc_G(G);
            cc_G.setEps(EpsRel);
            cc_G.setTol(EpsAbs);
            cc_G.setLim(Limit);
            cc_G.setRec(Recurrence);
            cc_G.setSubdiv(NNest);
            cc_G.setStack(NStack);
            cc_G.setRange(L);
            cc_G.setVerbose(Verbose, vName, Log);
            
            // integrate
            return cc_G.integrate_ff(-1., 1., n);    // (x1,+∞)->(-1,1)
        }

        // both bounds are finite
        return integrate_ff(x1, x2, n);
    }
    
    /**
     * @brief Finite integration.
     * 
     * Clenshaw-Curtis quadrature for finite interval (a,b).
     * 
     * @param x1 Left bound (allowed infinite).
     * @param x2 Right bound (allowed infinite).
     * @param n On output, evaluations needed for a converged result.
     */
    FType integrate_ff (double x1, double x2, int * n = nullptr) const
    {
        // check interval bounds
        if (x1 == x2)
        {
            if (n != nullptr)
                *n = 0;
            return FType(0);
        }
        if (x1 > x2)
        {
            return -integrate_ff(x2, x1, n);
        }
        
        // scaled F
        auto f = [&](double x) -> FType
        {
            return F(x1 + 0.5 * (1.0 + x) * (x2 - x1));
        };
        
        // Chebyshev coefficients
        std::vector<FType> coefs;
        
        // evaluated function f and its transformation
        std::vector<Complex> fvals_prev, fvals, ftraf;
        
        // integral approximation
        FType sum = 0, sum_prev = special::constant::Nan;
        
        // get nesting limit
        int maxN = gsl_sf_pow_int(2,NNest);
        
        // convergence loop
        for (int N = 4; N <= maxN /*or not Recurrence*/; N *= 2)
        {
            if (n != nullptr)
                *n = N;
            
            // reserve memory
            fvals.resize(2*N + 1);
            ftraf.resize(2*N);
            
            // create plan (thread-safe)
            fftw_plan plan;
#ifdef _OPENMP
            # pragma omp critical
#endif
            plan = fftw_plan_dft_1d
            (
                2*N,
                reinterpret_cast<fftw_complex*>(&fvals[0]),
                reinterpret_cast<fftw_complex*>(&ftraf[0]),
                FFTW_BACKWARD,
                0
            );
            
            // is this the first iteration?
            if (coefs.empty())
            {
                // evaluate f everywhere
                double pi_over_N = special::constant::pi / N;
                for (int k = 0; k < N; k++)
                {
                    fvals[k] = fvals[2*N-k] = f(std::cos(k * pi_over_N));
                    
                    if (not std::isfinite(std::abs(fvals[k])))
                        HexException("%s \"%g\" when evaluating function at %g", vName.c_str(), std::abs(fvals[k]), std::cos(k * pi_over_N));
                }
                fvals[N] = f(-1);
                
                if (not std::isfinite(std::abs(fvals[N])))
                    HexException("%s \"%g\" when evaluating function at -1.", vName.c_str(), std::abs(fvals[N]));
            }
            else
            {
                // evaluate just the new half, recycle older evaluations
                double pi_over_N = special::constant::pi / N;
                for (int k = 0; k < N; k++)
                {
                    fvals[k] = fvals[2*N-k] = (k % 2 == 0) ? fvals_prev[k/2] : f(std::cos(k * pi_over_N));
                    
                    if (not std::isfinite(std::abs(fvals[k])))
                        HexException("%s \"%g\" when evaluating function.", vName.c_str(), std::abs(fvals[k]));
                }
                fvals[N] = fvals_prev[N/2];
            }
            
            // append element
            coefs.resize(N + 1);

            // compute coefficients using FFT/DCT-I
            fftw_execute(plan);
            
            // delete plan (thread-safe)
#ifdef _OPENMP
            # pragma omp critical
#endif
            fftw_destroy_plan(plan);
            
            // create type-correct pointer
            FType const * ftraf_ptr = reinterpret_cast<FType*>(&ftraf[0]);
            
            // copy result
            if (typeid(FType) == typeid(Complex))
            {
                // copy whole complex numbers
                for (int i = 0; i <= N; i++)
                    coefs[i] = 0.5 * (*(ftraf_ptr + i));
            }
            else if (typeid(FType) == typeid(double))
            {
                // copy just real parts
                for (int i = 0; i <= N; i++)
                    coefs[i] = 0.5 * (*(ftraf_ptr + 2*i));
            }
            else
            {
                HexException("%s Can't handle datatype \"%s\".", vName.c_str(), typeid(FType).name());
            }
            
            // sum the quadrature rule
            sum = 0.5 * (coefs[0] - coefs[N] / (N*N - 1.));
            for (int twok = 2; twok < N; twok += 2)
                sum -= coefs[twok] / (twok*twok - 1.);
            
            // echo debug information
            if (Verbose)
                Log << vName << " N = " << N << ", Sum = " << FType(2.*(x2-x1)/N)*sum << "\n";
            
            // check for convergence
            if (std::abs(sum - FType(2.) * sum_prev) <= std::max(EpsRel*std::abs(sum), EpsAbs))
            {
                if (Verbose)
                    Log << vName << " Convergence for N = " << N << ", sum = " << FType(2. * (x2 - x1) / N) * sum << "\n";
                
                return FType(2. * (x2 - x1) / N) * sum;
            }
            else if ( std::isfinite(std::abs(sum)) 
                  and std::isfinite(std::abs(sum_prev)) 
                  and std::max(std::abs(sum), std::abs(sum_prev)) <= EpsAbs * std::abs(x2-x1) )
            {
                if (Verbose)
                    Log << vName << " EpsAbs limit matched, " << EpsAbs << " on (" << x1 << "," << x2 << ").\n";
                
                return FType(0.);
            }
            else
            {
                /* do nothing */
            }

            // save function evaluations and sum
            fvals_prev = std::move(fvals);
            sum_prev = sum;
        }
        
        // At this point the subdivision loop has been exited.
        // Now we either return, recurr or throw, depending on the settings.
        
        if (not Recurrence)
        {
            if (Throw)
            {
                HexException("%s Insufficient evaluation limit %d", vName.c_str(), maxN);
            }
            else
            {
                Log << vName << " WARNING: Insufficient evaluation limit " << maxN << ".\n";
                return FType(2. * (x2 - x1) / maxN) * sum;
            }
        }
        
        //
        // no convergence? -> bisect
        //
        
        // cancel bisection if interval too tiny
        if (std::abs(x2-x1) < EpsAbs)
        {
            if (Verbose)
                Log << vName << " Interval smaller than " << EpsAbs << "\n";
            return 0;
        }
        
        // cancel bisection if stack full
        if (NStack == 0)
        {
            if (Verbose)
                Log << vName << " Bisection inhibited due to internal stack limit.\n";
            return FType(2. * (x2 - x1) / maxN) * sum;
        }
        
        if (Verbose)
        {
            Log << vName << " Bisecting to ("
                << x1 << "," << (x2+x1)/2 << ") and ("
                << (x2+x1)/2 << "," << x2 << ")\n";
        }
        
        // the actual bisection - also the attributes of the class need to be modified
        
        int n1, n2;
        double this_EpsAbs = EpsAbs;
        
        NStack--;
        EpsAbs = 0.5 * (2. * (x2 - x1) / maxN) * std::abs(sum) * EpsRel;
        
        FType i1 = integrate_ff(x1, (x2+x1)/2, &n1);
        FType i2 = integrate_ff((x2+x1)/2, x2, &n2);
        
        NStack++;
        EpsAbs = this_EpsAbs;
        
        if (n != nullptr) *n = n1 + n2;
        return i1 + i2;
    }
    
    /**
     * @brief Improper integration.
     * 
     * Clenshaw-Curtis quadrature for infinite-infinite interval (-∞,+∞).
     * 
     * @param n On output, evaluations needed for converged result.
     */
    FType integrate_ii (int * n = nullptr) const
    {
        // function values, new and previous
        std::vector<FType> fvals, fvals_prev;

        // weights, new and previous
        std::vector<double> weights, weights_prev;

        // previous integral
        FType sum_prev = special::constant::Nan;

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
                    double x = i * special::constant::pi / N;
                    fvals[i] = F(L / std::tan(x));
                    if (not std::isfinite(std::abs(fvals[i])))
                        fvals[i] = 0;
                    weights[i] = 1. / gsl_sf_pow_int(std::sin(x), 2);
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
                        double x = i * special::constant::pi / N;
                        fvals[i] = F(L / std::tan(x));
                        if (not std::isfinite(std::abs(fvals[i])))
                            fvals[i] = 0;
                        weights[i] = 1. / gsl_sf_pow_int(std::sin(x), 2);
                    }
                }
            }

            // evaluate the integral
            FType sum = 0.;
            for (int i = 1; i <= N - 1; i++)
                sum += weights[i] * fvals[i];
            if (std::abs(sum - FType(2.) * sum_prev) < EpsRel * std::abs(sum))
            {
                if (n != nullptr) *n = N;
                return FType(L * special::constant::pi / N) * sum;
            }

            // save precomputed values
            sum_prev = sum;
            fvals_prev = fvals;
            weights_prev = weights;
        }
    }
    
private:
    
    /// Function to integrate
    Functor const & F;
    
    /// Relative precision
    double EpsRel;
    
    /// Relative precision
    mutable double EpsAbs;
    
    /**
     * @brief Whether to use @ref lim for evaluating integration boundaries.
     * 
     * Whether to use limit (limit = true) when evaluating F(x1) and
     * F(x2) for improper arguments x1=-∞ and/or x2=+∞. Otherwise 
     * (limit = false) the functor will be simply evaluated at
     * possibly infinite boundaries and it has to cope itself with the
     * input.
     */
    bool Limit;
    
    /**
     * @brief Enable recurrent subdivision.
     * 
     * Whether to divide-and-conquer if current level doesn't converge
     * after NLevel evaluations.
     */
    bool Recurrence;
    
    /**
     * @brief Breadth limit for bisection.
     * 
     * Sets the maximal evaluations count for a single bisection level.
     * When the limit is reached, a new bisection will be done. Thus, the
     * integration rule will be nested NLevel-times before doing so.
     */
    int NNest;
    
    /**
     * @brief Depth limit for bisection.
     * 
     * Sets the subdivision level count. The initial integration interval
     * will be subdivided into 2^NStack pieces if totally non-convergent.
     */
    mutable int NStack;
    
    /// Compactification parameter
    double L;
    
    /// Display debugging information.
    bool Verbose;
    
    /// Debuggin information identification.
    std::string vName;
    
    /// Debugging stream (defaults to 'std::cout').
    mutable std::ostream Log;
    
    /// Throw on non-critical errors.
    bool Throw;
};

#endif
