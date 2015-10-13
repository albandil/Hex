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

#ifndef HEX_COMPACT
#define HEX_COMPACT

#include <cassert>
#include <cmath>
#include <iostream>

#include "hex-misc.h"
#include "hex-special.h"

/**
 * @brief Compute limit.
 * 
 * Limit
 * @f[
 *     \lim_{t \rightarrow x} F(t)
 * @f]
 *
 * The bisection count will be written into *n if non-null.
 */
template <class Functor, typename FType> FType lim (Functor F, double x, int * n = nullptr)
{
    // initial position
    double x0 = 0;
    if (x > 0.) x0 = std::isfinite(x) ? 0.5 * x :  1.;
    if (x < 0.) x0 = std::isfinite(x) ? 0.5 * x : -1.;

    // function value
    FType f0 = F(x0), f = special::constant::Nan;

    // main loop
    int i;
    for (i = 0; n == nullptr or i <= *n; i++)
    {
        // advance x0
        if (x == special::constant::Inf)
            x0 *= 2.;
        else if (x == -special::constant::Inf)
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
 * @brief Compactification.
 * 
 * An auxiliary interface class. All other compactification classes
 * are derived from this class.
 */
template <typename FType> class ICompactification
{
public:
    virtual ~ICompactification() {}
    virtual double scale (double x) const = 0;
    virtual double unscale (double t) const = 0;
    virtual double Jacobian (double t) const = 0;
    virtual FType operator() (double t) const = 0;
};

/**
 * @brief Compactification of a function from finite interval.
 * 
 * Transform arbitrary function so that its definition range will
 * be @f$ t \in [-1,1] @f$, for original interval @f$ x \in [a,b] @f$. The
 * following formula is used:
 *
 * @f[
 *     x = \frac{a+b}{2} + t \frac{b-a}{2}
 *     \Leftrightarrow
 *     t = \frac{x - (a+b)/2}{(b-a)/2} \ .
 * @f]
 */
template <class Functor, typename FType> class CompactificationF : public ICompactification<FType>
{
public:
    /**
     * @brief Constructor of the class.
     * The parameters "a" and "b" specify original definition interval [a,b].
     * It must be a < b.
     */
    CompactificationF (
        Functor f,
        double a,
        double b
    ) : F(f), A(a), B(b), M(0.5*(b+a)), D(0.5*(b-a))
    {
        if (not std::isfinite(a) or not std::isfinite(b))
            HexException("[CompactificationF] Interval has to be finite!");
    }

    /// Scale value from the original interval [a,b] into compactified interval [-1,1].
    double scale (double x) const
    {
        assert(A <= x and x <= B);
        return (x - M) / D;
    }

    /// Unscale value from the compactified interval [-1,1] into the original interval [a,b].
    double unscale (double t) const
    {
        assert(std::abs(t) <= 1.);
        return M + t * D;
    }

    /// Evaluate Jacobian of the transformation.
    double Jacobian(double t) const
    {
        assert(std::abs(t) <= 1.);
        return D;
    }

    /**
     * @brief Evaluate the compactified function.
     * @param t Value from the compactified interval [-1,1].
     */
    FType operator() (double t) const
    {
        return F(M + t * D);
    }

private:
    Functor F;
    double A, B, M, D;
};

/**
 * @brief Compactification of a function from (-∞,b] to (-1,1].
 * 
 * Transform arbitrary function so that its definition range will
 * be @f$ t \in [-1,1] @f$, for original interval @f$ x \in [a,b] @f$. The
 * following formulas are used:
 *
 * @f[
 *     x = b - L \frac{1+t}{1-t}
 *     \Leftrightarrow
 *     t = \frac{x - b + L}{x - b + L} \ .
 * @f]
 *   Note that in this case the monotonic character is reverse.
 */
template <class Functor, typename FType> class CompactificationL : public ICompactification<FType>
{
public:

    /**
     * Constructor of the class
     * The parameter "b" specifies original definition interval (-∞,b].
     */
    CompactificationL (
        Functor f,
        double b = 0.,
        bool limit = true,
        double L = 1.0
    ) : F(f), B(b), L(L), Limit(limit) {}

    /// Scale value from the original interval [a,b] into compactified interval [-1,1].
    double scale (double x) const
    {
        assert (x <= B);
        return std::isfinite(x) ? (x - B + L) / (x - B - L) : 1.;
    }

    /// Unscale value from the compactified interval [-1,1] into the original interval [a,b].
    double unscale (double t) const
    {
        assert(std::abs(t) <= 1.);
        return (t == 1.) ? special::constant::Inf : B - L * (1. + t) / (1. - t);
    }

    /// Evaluate Jacobian of the transformation.
    double Jacobian(double t) const
    {
        assert(std::abs(t) <= 1.);
        return -2 * L / ((1. - t) * (1. - t));
    }

    /**
     * Evaluate the compactified function.
     * @param t Value from the compactified interval [-1,1].
     */
    FType operator() (double t) const
    {
        if (t == 1.)
        {
            if (Limit)
                return lim<decltype(F),FType>(F, special::constant::Inf);
            else
                return F(special::constant::Inf);
        }
        else
        {
            return F(B - L * (1. + t) / (1. - t));
        }
    }

private:
    Functor F;
    double B, L;
    bool Limit;
};

/**
 * @brief Compactification of a function from [a,+∞) to [-1,1).
 * 
 * Transform arbitrary function so that its definition range will
 * be @f$ t \in [-1,1] @f$, for original interval @f$ x \in [a,b] @f$. The
 * following formulas are used:
 *
 * @f[
 *     x = a + L \frac{1+t}{1-t}
 *     \Leftrightarrow
 *     t = \frac{x - a - L}{x - a - L} \ .
 * @f]
 *
 * - Compactification of a double-unbounded interval is not implemented yet.
 */
template <class Functor, typename FType> class CompactificationR : public ICompactification<FType>
{
public:

    /**
     * Constructor of the class
     * The parameter "a" specifies original definition interval @f$ [a,+\infty) @f$.
     */
    CompactificationR (
        Functor f,
        double a = 0.,
        bool limit = true,
        double L = 1.0
    ) : F(f), A(a), L(L), Limit(limit) {}

    /// Scale value from the original interval [a,b] into compactified interval [-1,1].
    double scale (double x) const
    {
        assert (x >= A);
        return std::isfinite(x) ? (x - A - L) / (x - A + L) : 1.;
    }

    /// Unscale value from the compactified interval [-1,1] into the original interval [a,b].
    double unscale (double t) const
    {
        assert(std::abs(t) <= 1.);
        return (t == 1.) ? special::constant::Inf : A + L * (1. + t) / (1. - t);
    }

    /// Evaluate Jacobian of the transformation.
    double Jacobian(double t) const
    {
        assert(std::abs(t) <= 1.);
        return 2 * L / ((1. - t) * (1. - t));
    }

    /**
     * Evaluate the compactified function.
     * @param t Value from the compactified interval [-1,1].
     */
    FType operator() (double t) const
    {
        if (t == 1.)
        {
            if (Limit)
                return lim<decltype(F),FType>(F,special::constant::Inf);
            else
                return F(special::constant::Inf);
        }
        else
        {
            return F(A + L * (1. + t) / (1. - t));
        }
    }

private:
    Functor F;
    double A, L;
    bool Limit;
};

/**
 * @brief Compactification multiplied by its jacobian.
 * 
 * A wrapper around the compactification classes which returns function
 * value multiplied by the Jacobian of the compactification transformation.
 * It is meant for the use in integration, so that one can just call its
 * operator() interface instead of operator()*jacobian
 */
template <class Functor, typename FType> class CompactIntegrand
{
public:

    /**
     * Constructor of the class
     * The parameters "a" and "b" specify original definition interval [a,b].
     * It must be a < b.
     */
    CompactIntegrand (
        Functor f,
        double a = 0.,
        double b = special::constant::Inf,
        bool limit = true,
        double L = 1.0
    ) : Compactification(nullptr) {
        if (std::isfinite(a) and std::isfinite(b))
            Compactification = new CompactificationF<Functor,FType> (f, a, b);
        else if (std::isfinite(a) and not std::isfinite(b))
            Compactification = new CompactificationR<Functor,FType> (f, a, limit, L);
        else if (not std::isfinite(a) and std::isfinite(b))
            Compactification = new CompactificationL<Functor,FType> (f, b, limit, L);
        else
            HexException("[CompactIntegrand] Compactification of (-∞,∞) interval is not implemeted.");
    }

    ~CompactIntegrand ()
    {
        delete Compactification;
    }

    /**
     * Get the function multiplied by the Jacobian of the transform for the
     * purpose of integrating the original function.
     *
     * @note Whenever the function value is zero, the Jacobian is not evaluated
     *       and a clean zero is returned.
     */
    inline FType operator() (double t) const
    {
        // evaluate function
        FType ft = Compactification->operator()(t);

        if (ft == FType(0.))
        {
            // return zero
            return 0.;
        }
        else
        {
            // multiply by the Jacobian
            return ft * Compactification->Jacobian(t);
        }
    }

    /// Scale value from the original interval [a,b] into compactified interval [-1,1].
    inline double scale (double x) const
    {
        return Compactification->scale(x);
    }

    /// Unscale value from the compactified interval [-1,1] into the original interval [a,b].
    inline double unscale (double t) const
    {
        return Compactification->unscale(t);
    }

private:

    ICompactification<FType> * Compactification;
};

#endif
