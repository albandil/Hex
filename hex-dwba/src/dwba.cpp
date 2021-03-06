//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#include "hex-arrays.h"
#include "hex-hydrogen.h"
#include "hex-gausskronrod.h"
#include "hex-special.h"

#include "potential.h"
#include "dwba.h"
#include "wave_distort.h"

using special::constant::pi;

namespace DWBA1
{

Complex computeDirect1e (DistortingPotential const& U, int l, double k)
{
    // get distorted wave
    DistortedWave chi_kl(k,l,U);

    // set up the integrands
    auto integrand1 = [ & ](double r) -> double
    {
        return chi_kl(r) * U(r) * special::ric_j(l,k*r);
    };
    auto integrand2 = [ & ](double r) -> double
    {
        return chi_kl(r) * U(r) * chi_kl(r);
    };

    // integrate
    GaussKronrod<decltype(integrand1)> Q1(integrand1);
    GaussKronrod<decltype(integrand2)> Q2(integrand2);
    Q1.integrate(0., special::constant::Inf);
    Q2.integrate(0., special::constant::Inf);

    // return the result
    return std::pow(4*pi,2) * std::sqrt((2*l+1)/(4*pi)) / (k*k) * 
        (Q1.result() - Q2.result() * chi_kl.getPhasef()) * chi_kl.getPhasef();
}

Complex computeExchange1e
(
    DistortingPotential const& U,
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf
)
{
    // get distorted waves
    DistortedWave chif(kf,Li,U);
    DistortedWave chii(ki,Lf,U);

    // set up the integrands
    auto integrand1 = [ & ](double r) -> double
    {
        return chif(r) * U(r) * Hydrogen::P(Ni,Li,r);
    };
    auto integrand2 = [ & ](double r) -> double
    {
        return Hydrogen::P(Nf,Lf,r) * chii(r);
    };

    // integrate
    GaussKronrod<decltype(integrand1)> Q1(integrand1);
    GaussKronrod<decltype(integrand2)> Q2(integrand2);
    Q1.integrate(0., special::constant::Inf);
    Q2.integrate(0., special::constant::Inf);

    // compute phase factor
    double phase = chif.getPhase() + chii.getPhase();
    Complex phasef(cos(phase),sin(phase));

    // return the result
    return -pow(4*pi,2) * phasef * pow(Complex(0.,1.), Lf-Li) * 
        sqrt((2*Lf+1)/(4*pi)) * Q1.result() * Q2.result() / (ki*kf);
}

Complex computeDirect2e
(
    const DistortingPotential& U, int lambda,
    int Nf, int Lf, double kf, int lf,
    int Ni, int Li, double ki, int li
)
{
    // get distorted waves
    DistortedWave chif(kf,lf,U);
    DistortedWave chii(ki,li,U);

    // set up the outer integrand
    auto outer_integrand = [ & ](double r2) -> double
    {
        // set up the inner integrands
        auto inner_integrand_1 = [ & ](double r1) -> double
        {
            double psi_f = Hydrogen::P(Nf,Lf,r1);
            double psi_i = Hydrogen::P(Ni,Li,r1);
            double multipole;
            if (lambda == 0)
                multipole = 1/r1 - 1/r2;
            else
                multipole = pow(r2/r1,lambda)/r1;
            return psi_f * multipole * psi_i;
        };
        auto inner_integrand_2 = [ & ](double r1) -> double
        {
            double psi_f = Hydrogen::P(Nf,Lf,r1);
            double psi_i = Hydrogen::P(Ni,Li,r1);
            double multipole;
            if (lambda == 0)
                return 0;
            else
                multipole = pow(r1/r2,lambda)/r2;
            return psi_f * multipole * psi_i;
        };

        // integrate
        GaussKronrod<decltype(inner_integrand_1)> Q1(inner_integrand_1);
        GaussKronrod<decltype(inner_integrand_2)> Q2(inner_integrand_2);
        Q1.integrate(r2, special::constant::Inf);
        Q2.integrate(0, r2);

        // evaluate distorted waves
        double chi_i = chii(r2);
        double chi_f = chif(r2);

        // return the result
        return chi_i * chi_f * (Q1.result() + Q2.result());
    };

    // integrate
    GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
    Q.integrate(0., special::constant::Inf);

    // compute phase factor
    double phase = chii.getPhase() + chif.getPhase();
    Complex phasef(cos(phase),sin(phase));

    // return the result
    return pow(4*pi,2) * phasef * pow(Complex(0.,1.),li-lf) * 
        sqrt((2*li+1)/(4*pi)) * Q.result() / (ki*kf);
}

Complex computeExchange2e
(
    const DistortingPotential& U, int lambda,
    int Nf, int Lf, double kf, int lf,
    int Ni, int Li, double ki, int li
)
{
    // get distorted waves
    DistortedWave chif(kf,lf,U);
    DistortedWave chii(ki,li,U);

    // set up the outer integrand
    auto outer_integrand = [ & ](double r2) -> double {

        // set up the inner integrands
        auto inner_integrand_1 = [ & ](double r1) -> double {
            double psi_f = Hydrogen::P(Nf,Lf,r1);
            double chi_i = chii(r1);
            double multipole;
            if (lambda == 0)
                return 0.;
            else
                multipole = pow(r2/r1,lambda)/r1;
            return psi_f * multipole * chi_i;
        };
        auto inner_integrand_2 = [ & ](double r1) -> double {
            double psi_f = Hydrogen::P(Nf,Lf,r1);
            double chi_i = chii(r1);
            double multipole;
            if (lambda == 0)
                multipole = 1/r2 - 1/r1;
            else
                multipole = pow(r1/r2,lambda)/r2;
            return psi_f * multipole * chi_i;
        };

        // integrate
        GaussKronrod<decltype(inner_integrand_1)> Q1(inner_integrand_1);
        GaussKronrod<decltype(inner_integrand_2)> Q2(inner_integrand_2);
        Q1.integrate(r2, special::constant::Inf);
        Q2.integrate(0, r2);

        // evaluate distorted waves
        double psi_i = Hydrogen::P(Ni,Li,r2);
        double chi_f = chif(r2);

        // return the result
        return psi_i * chi_f * (Q1.result() + Q2.result());
    };

    // integrate
    GaussKronrod<decltype(outer_integrand)> Q(outer_integrand);
    Q.integrate(0., special::constant::Inf);

    // compute phase factor
    double phase = chii.getPhase() + chif.getPhase();
    Complex phasef(cos(phase),sin(phase));

    // return the result
    return pow(4*pi,2) * phasef * pow(Complex(0.,1.),li-lf) * 
        sqrt((2*li+1)/(4*pi)) * Q.result() / (ki*kf);
}

} // end of namespace DWBA1

void dwba
(
    int Ni, int Li, double ki,
    int Nf, int Lf, double kf,
    int L,
    cArrays & Tdir, cArrays & Texc,
    double rmax, bool direct, bool exchange
)
{
    // m⟶m' transition count
    int MM = (2*Li+1)*(2*Lf+1);

    // allocate memory
    Tdir.resize(MM);
    Texc.resize(MM);

    // initial and final atomic state
    HydrogenFunction psii(Ni, Li);
    HydrogenFunction psif(Nf, Lf);

    // distorting potentials
    DistortingPotential Uf(Nf,rmax);
    DistortingPotential Ui(Ni,rmax);

    std::cout << "L = " << L << std::endl;

    for (int lf = std::abs(Lf - L); lf <= Lf + L; lf++)
    {
        // add new T-matrix for this outgoing partial wave
        for (cArray & T : Tdir)
            T.push_back(0.);
        for (cArray & T : Texc)
            T.push_back(0.);

        DistortedWave chif(kf,lf,Ui);

        // direct 1e
        if (Ni == Nf and Li == Lf and direct)
        {
            Complex tmat = DWBA1::computeDirect1e(Uf,lf,ki);

            std::cout << "\tlf = " << lf << ", li = " << lf << " : direct 1e = " << tmat << std::endl;

            for (int Mi = -Li; Mi <= Li; Mi++)
            for (int Mf = -Lf; Mf <= Lf; Mf++)
                Tdir[(Mi+Li)*(2*Lf+1)+Mf+Lf][lf-std::abs(Lf - L)] += tmat;
        }

        for (int li = std::abs(Li - L); li <= Li + L; li++)
        {
            // conserve angular momentum
            if (li < lf - Li - Lf)
                continue;

            // conserve angular momentum
            if (li > lf + Li + Lf)
                break;

            // conserve parity
            if ((li + Li) % 2 != (lf + Lf) % 2)
                continue;

            cArray DD_lf_li(MM), DE_lf_li(MM), ED_lf_li(MM), EE_lf_li(MM);
            DistortedWave chii(ki,li,Ui);

            // exchange 1e
            if (Li == lf and Lf == li and exchange)
            {
                Complex tmat = DWBA1::computeExchange1e(Uf, Ni, Li, ki, Nf, Lf, kf);

                std::cout << "\tlf = " << lf << ", li = " << li << " : exchange 1e = " << tmat << std::endl;

                for (int Mi = -Li; Mi <= Li; Mi++)
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                        Texc[(Mi+Li)*(2*Lf+1)+Mf+Lf][lf-std::abs(Lf - L)] += (Mf == 0) ? tmat : 0.;
            }

            // direct 2e
            for (int lambda = std::max(abs(Li-Lf),abs(li-lf)); direct and lambda <= std::min(Li+Lf,li+lf); lambda++)
            {
                Complex tmat = DWBA1::computeDirect2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);

                std::cout << "\tlf = " << lf << ", li = " << li << ", λ = " << lambda << " : direct 2e = " << tmat << std::endl;

                for (int Mi = -Li; Mi <= Li; Mi++)
                {
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                    {
                        double ang = special::ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi)
                                   * special::ClebschGordan(Li,Mi,li,0,L,Mi)
                                   * special::computef(lambda,Lf,lf,Li,li,L);

                        if (not std::isfinite(ang))
                            throw exception ("Gaunt failure!\n");

                        Tdir[(Mi+Li)*(2*Lf+1)+Mf+Lf][lf-std::abs(Lf - L)] += tmat * ang;
                    }
                }
            }

            // exchange 2e
            for (int lambda = std::max(abs(Li-lf),abs(li-Lf)); exchange and lambda <= std::min(Li+lf,li+Lf); lambda++)
            {
                Complex tmat = DWBA1::computeExchange2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);

                std::cout << "\tlf = " << lf << ", li = " << li << ", λ = " << lambda << " : exchange 2e = " << tmat << std::endl;

                for (int Mi = -Li; Mi <= Li; Mi++)
                {
                    for (int Mf = -Lf; Mf <= Lf; Mf++)
                    {
                        double ang = special::ClebschGordan(Lf,Mf,lf,Mi-Mf,L,Mi)
                                   * special::ClebschGordan(Li,Mi,li,0,L,Mi)
                                   * special::computef(lambda,lf,Lf,Li,li,L);

                        if (not std::isfinite(ang))
                            throw exception ("Gaunt failure!\n");

                        Texc[(Mi+Li)*(2*Lf+1)+Mf+Lf][lf-std::abs(Lf - L)] += tmat * ang;
                    }
                }
            }
        } /* for li */
    } /* for lf */
}
