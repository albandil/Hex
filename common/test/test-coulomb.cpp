#include <cmath>
#include <cstdlib>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

#include "hex-special.h"

// Test of the regular Coulomb wave function of an attractive field (Z = 1).
//
// The point of interest is the region of very small wavenumbers. There the Sommerfeld
// parameter η = -Z/k is enormous, the continued fractions used by GSL stop converging
// and gsl_sf_coulomb_wave_FG_e returns GSL_ERUNAWAY (or, for large angular momenta,
// GSL_ELOSS). The function itself is perfectly well behaved there, so special::coul_F has
// to supply the value on its own; it does so from the zero-energy limit
//
//     F_l(-Z/k, kr)  ⟶  sqrt(π k r) J_{2l+1}(sqrt(8 Z r)) .
//
// The reference values below come from mpmath 1.3.0 (mp.dps = 30, mp.coulombf and its
// derivative with respect to ρ = kr). The last six of them lie in the region GSL refuses.
// An earlier version of the code answered the runaway ones with the value at r = 0, which
// is not merely imprecise but unrelated to the function at r > 0, F being oscillatory
// however small the energy is, and the loss-of-accuracy ones with the uniform
// approximation, which is off by per cent this far from the classical turning point.
// That is what test 1 guards against.
//
// A loss of accuracy is also reported in a second, unrelated place: below the classical
// turning point of a large angular momentum, where F is exponentially small. GSL cannot
// be trusted there either -- over a sweep of that region its value came out non-finite in
// a third of the cases -- so coul_F answers from the uniform approximation, which holds to
// a couple of per cent throughout. Test 6 covers that region; note that the zero-energy
// limit is of no use in it, k**2 r being no longer small.
//
// The whole thing takes a fraction of a second, so like the sparse grid test it can be
// run after every build.

struct Reference
{
    int l;
    double k, r, F, Fp;
    char const * comment;
};

const Reference references [] =
{
    { 0, 1.0,        5.0,      0.90941019611717469,   0.1775696135827057,   "GSL is comfortable here" },
    { 3, 0.5,        20.0,     0.016803811739245932, -1.0634110290964349,   "GSL is comfortable here" },
    { 0, 1.0e-4,     43.4178, -0.017908923476996489,  24.849798134926704,   "just above the boundary; GSL still converges" },
    { 0, 1.73555e-5, 43.4178, -0.0074608491620017126, 59.649090383740144,   "the case reported from hex-ecs extraction" },
    { 1, 1.73555e-5, 43.4178,  0.0062537708835755501,-77.849402654114713,   "the same, higher partial wave" },
    { 3, 1.0e-5,     100.0,    0.00076942265745799184, 116.77408187111051,  "deep in the zero-energy regime" },
    { 6, 1.0e-6,     1000.0,   0.004743345634227142, -12.89736455451356,    "large radius, high partial wave" },
    { 0, 1.0e-6,     3000.0,  -0.0012261230442263072,-157.63713836779468,   "the extreme corner of the acceptance region" },
    { 6, 3.84434e-6, 10.0,      0.00011200257995791101, 15.745233053815797,  "small k announced by GSL as a loss of accuracy" }
};

// The zero-energy limit neglects k² next to 2Z/r, so its relative error grows with k²r.
// GSL abandons the evaluation while k²r ~ 1e-8, where the limit holds to some eight
// digits; the reference at r = 3000 sits at k²r = 3e-9 and is the least accurate of the
// set. The bound is loose enough for that and still tight enough to notice any change
// of the formula or of the branch that selects it.
const double tolerance = 1e-6;

// Points below the classical turning point of a large angular momentum, ordered by
// decreasing rho/rho_t. GSL reports a loss of accuracy at every one of them and coul_F
// falls back on the uniform approximation, whose error over this region was measured
// against mpmath at 0.4 to 2.1 per cent in F and up to 4.3 per cent in F'; the tolerance
// below leaves room for that and still catches a return to GSL's value, to the
// zero-energy limit, or to anything else that is wrong by more than a small factor.
const Reference barriers [] =
{
    { 15, 0.0650239, 25.0365, 9.26781313336667025e-08, 8.12040353479437451e-07, "the case reported from hex-ecs extraction; rho/rho_t = 0.25" },
    { 30, 0.05,      60.0,    4.21581649982650681e-18, 4.05398881463497470e-17, "rho/rho_t = 0.18; GSL's own value is 70 % off here" }
};

// See the comment above the table.
const double barrier_tolerance = 5e-2;

bool check (char const * what, double result, double expected, double tol)
{
    double error = std::abs(result - expected) / std::max(std::abs(expected), 1e-300);

    std::cout << "\t" << what << " = " << result << " (expected " << expected << ")"
              << ", rel. error " << error << (error < tol ? "" : "  <-- above the tolerance")
              << std::endl;

    return error < tol;
}

int main (void)
{
    // GSL must not abort the program on the errors that are handled internally
    gsl_set_error_handler_off();

    bool passed = true;

    // 1. values against the high-precision references

    std::cout << "Test 1: values of coul_F" << std::endl;

    for (Reference const & ref : references)
    {
        std::cout << "  l = " << ref.l << ", k = " << ref.k << ", r = " << ref.r
                  << " (" << ref.comment << ")" << std::endl;

        double F, Fp;
        int err = special::coul_F(1, ref.l, ref.k, ref.r, F, Fp);

        if (err != GSL_SUCCESS)
        {
            std::cout << "\tcoul_F failed with " << gsl_strerror(err) << std::endl;
            passed = false;
            continue;
        }

        passed &= check("F ", F,  ref.F,  tolerance);
        passed &= check("F'", Fp, ref.Fp, tolerance);
    }

    // 2. scaling with the wavenumber
    //    At vanishing energy the radial shape no longer depends on k and the whole
    //    dependence is the sqrt(k) of the normalization. This holds no matter what
    //    the value of the function is, so it also catches an error common to all of
    //    the references above.

    std::cout << "Test 2: F is proportional to sqrt(k) at vanishing energy" << std::endl;

    for (int l = 0; l <= 4; l += 2)
    {
        double F1, Fp1, F2, Fp2;
        special::coul_F(1, l, 1.0e-5, 137.0, F1, Fp1);
        special::coul_F(1, l, 4.0e-6, 137.0, F2, Fp2);

        char what [] = "F(k)/F(k/2.5), l = 0";
        what[sizeof(what) - 2] = '0' + l;

        passed &= check(what, F1 / F2, std::sqrt(2.5), tolerance);
    }

    // 3. the derivative against a finite difference of the value
    //    Both are supplied by the same branch of coul_F, so this does not repeat the
    //    reference check; it verifies that the derivative is taken with respect to
    //    ρ = kr, which is what the callers assume, and not with respect to r.

    std::cout << "Test 3: F' is the derivative with respect to rho = kr" << std::endl;

    for (int l = 0; l <= 2; l++)
    {
        const double k = 2.0e-5, r = 60.0, h = 1.0e-3;

        double Fm, Fc, Fp_, dummy, Fprho;
        special::coul_F(1, l, k, r - h, Fm,  dummy);
        special::coul_F(1, l, k, r,     Fc,  Fprho);
        special::coul_F(1, l, k, r + h, Fp_, dummy);

        // dF/dρ = (dF/dr) / k
        double difference = (Fp_ - Fm) / (2 * h * k);

        char what [] = "dF/drho, l = 0";
        what[sizeof(what) - 2] = '0' + l;

        // the central difference is second order accurate, hence the looser bound
        passed &= check(what, Fprho, difference, 1e-5);

        // silence the warning about the unused value
        (void)Fc;
    }

    // 4. the uniform approximation must refuse arguments it cannot handle
    //    It is built around the classical turning point, which an attractive field
    //    does not have for l = 0. It used to return NaN there, which coul_F would
    //    have passed on to the caller had GSL ever reported a loss of accuracy.

    std::cout << "Test 4: coul_F_michel refuses l = 0 in an attractive field" << std::endl;

    double Fm, Fpm;
    int errm = special::coul_F_michel(1, 0, 0.1, 5.0, Fm, Fpm);

    // NOTE: The refusal has to be the deliberate one. Without it the routine happens
    //       to report an underflow as well, that being what the Airy call makes of the
    //       NaN it is given, so accepting any failure here would accept the NaN too.
    std::cout << "\treturned " << gsl_strerror(errm)
              << (errm == GSL_EDOM ? "" : "  <-- expected a domain error") << std::endl;

    if (errm != GSL_EDOM)
        passed = false;

    // 5. the region GSL gives up on is the region that is handled internally
    //    This is not a property of the Coulomb function but of the two implementations
    //    meeting; it fails if a future GSL widens or narrows its range of convergence
    //    beyond the margin the acceptance criterion of coul_F leaves.

    std::cout << "Test 5: coul_F succeeds wherever GSL alone does not" << std::endl;

    int nrunaway = 0, nfailed = 0;

    for (double r = 10.0; r <= 3000.0; r *= 3.0)
    for (int l = 0; l <= 6; l += 3)
    for (double k = 1.0e-7; k < 1.0e-2; k *= 1.5)
    {
        gsl_sf_result f, g, fp, gp;
        double ef, eg;

        if (gsl_sf_coulomb_wave_FG_e(-1.0 / k, k * r, l, 0, &f, &fp, &g, &gp, &ef, &eg) == GSL_SUCCESS)
            continue;

        nrunaway++;

        double F, Fp;
        if (special::coul_F(1, l, k, r, F, Fp) != GSL_SUCCESS or not std::isfinite(F) or not std::isfinite(Fp))
            nfailed++;
    }

    std::cout << "\t" << nrunaway << " points refused by GSL, " << nfailed
              << " of them not recovered" << (nfailed == 0 ? "" : "  <-- expected none")
              << std::endl;

    if (nrunaway == 0)
    {
        std::cout << "\tGSL refused nothing at all  <-- the sweep no longer covers the case" << std::endl;
        passed = false;
    }

    passed &= (nfailed == 0);

    // 6. the classically forbidden region of a large angular momentum
    //    This is the other face of the loss of accuracy, and it has nothing to do with
    //    a small wavenumber: below the classical turning point F falls off
    //    exponentially and GSL abandons it while still handing back a number, one that
    //    is wrong by tens of per cent or not finite at all. The uniform approximation
    //    is what coul_F has to answer with there, so the test pins both the value and
    //    the branch that produced it.

    std::cout << "Test 6: the region below the classical turning point" << std::endl;

    for (Reference const & ref : barriers)
    {
        std::cout << "  l = " << ref.l << ", k = " << ref.k << ", r = " << ref.r
                  << " (" << ref.comment << ")" << std::endl;

        // the premise: GSL alone does not manage this point
        gsl_sf_result f, fpr, g, gpr;
        double expF, expG;
        int errgsl = gsl_sf_coulomb_wave_FG_e
        (
            -1.0 / ref.k, ref.k * ref.r, ref.l, 0,
            &f, &fpr, &g, &gpr, &expF, &expG
        );

        if (errgsl != GSL_ELOSS)
        {
            std::cout << "\tGSL returned " << gsl_strerror(errgsl)
                      << " instead of a loss of accuracy  <-- the point no longer tests"
                         " the fallback" << std::endl;
            passed = false;
        }

        double F, Fp;
        int err = special::coul_F(1, ref.l, ref.k, ref.r, F, Fp);

        if (err != GSL_SUCCESS)
        {
            std::cout << "\tcoul_F failed with " << gsl_strerror(err) << std::endl;
            passed = false;
            continue;
        }

        passed &= check("F ", F,  ref.F,  barrier_tolerance);
        passed &= check("F'", Fp, ref.Fp, barrier_tolerance);

        // The value has to be the one of the uniform approximation. Without this the
        // test would also accept the zero-energy limit at the first point, which is
        // only a factor of eight out and would slip through a per-cent tolerance at
        // some other radius.
        double Fu, Fpu;
        special::coul_F_michel(1, ref.l, ref.k, ref.r, Fu, Fpu);

        if (F != Fu or Fp != Fpu)
        {
            std::cout << "\tthe value did not come from the uniform approximation"
                         "  <-- expected coul_F_michel" << std::endl;
            passed = false;
        }
    }

    std::cout << (passed ? "All tests passed." : "Some tests failed.") << std::endl;

    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
