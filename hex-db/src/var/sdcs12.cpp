//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#include <map>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-chebyshev.h"
#include "hex-clenshawcurtis.h"
#include "hex-version.h"

// --------------------------------------------------------------------------------- //

#include "../quantities.h"
#include "../utils.h"

// --------------------------------------------------------------------------------- //

createNewScatteringQuantity(SingleDifferentialCrossSectionWrtRelativeAngle, "sdcs12")

// --------------------------------------------------------------------------------- //

std::string SingleDifferentialCrossSectionWrtRelativeAngle::description ()
{
    return "Single differential ionization cross section with respect to the relative scattering angle.";
}

std::vector<std::string> SingleDifferentialCrossSectionWrtRelativeAngle::dependencies ()
{
    return std::vector<std::string>
    {
        "ionf"
    };
}

std::vector<std::pair<std::string,std::string>> SingleDifferentialCrossSectionWrtRelativeAngle::params ()
{
    return std::vector<std::pair<std::string,std::string>>
    {
        {"ni", "Initial atomic principal quantum number."},
        {"li", "Initial atomic orbital quantum number."},
        {"mi", "Initial atomic magnetic quantum number."},
        {"S", "Total spin of atomic + projectile electron."},
        {"Ei", "Impact energy in Ry."},
        {"theta12", "Relative scattering angle, cos(theta12) = k1.k2 / |k1||k2|."}
    };
}

std::vector<std::string> SingleDifferentialCrossSectionWrtRelativeAngle::vparams ()
{
    return std::vector<std::string>
    {
        "theta12"
    };
}

// --------------------------------------------------------------------------------- //

bool SingleDifferentialCrossSectionWrtRelativeAngle::initialize (sqlitepp::session & db)
{
    return ScatteringQuantity::initialize(db);
}

bool SingleDifferentialCrossSectionWrtRelativeAngle::createTable ()
{
    return true;
}

bool SingleDifferentialCrossSectionWrtRelativeAngle::updateTable ()
{
    return ScatteringQuantity::updateTable();
}

// --------------------------------------------------------------------------------- //

bool SingleDifferentialCrossSectionWrtRelativeAngle::run (std::map<std::string,std::string> const & sdata)
{
    // manage units
    double efactor = change_units(Eunits, eUnit_Ry);
    double lfactor = change_units(lUnit_au, Lunits);
    double afactor = change_units(Aunits, aUnit_rad);

    // scattering event parameters
    int ni = Conv<int>(sdata, "ni", name());
    int li = Conv<int>(sdata, "li", name());
    int mi = Conv<int>(sdata, "mi", name());
    int  S = Conv<int>(sdata, "S",  name());
    double Ei = Conv<double>(sdata, "Ei", name()) * efactor;

    // maximal momentum
    double kmax = std::sqrt(Ei);

    // info structure
    typedef struct
    {
        int l1, l2;
        Chebyshev<double,Complex> cb;
    }
    Info;

    // parameters
    int l1, l2, L;
    std::string blob;
    std::map<int,std::vector<Info>> CB;
    rArray theta12;

    // get energy / energies
    try
    {
        // is there a single angle specified using command line ?
        theta12.push_back(Conv<double>(sdata, "theta12", name()));
    }
    catch (std::exception e)
    {
        // are there more angles specified using the STDIN ?
        theta12 = readStandardInput<double>();
    }

    // write header
    std::cout << logo("#") <<
        "# Single differential cross section (wrt relative angle) in " << unit_name(Lunits) << " for\n" <<
        "#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
        "#     S = " << S << ", Ei = " << Ei << " " << unit_name(Eunits) << "\n" <<
        "# ordered by relative angle in " << unit_name(Aunits) << "\n" <<
        "#\n";
    OutputTable table;
    table.setWidth(15);
    table.setAlignment(OutputTable::left);
    table.write("# theta12  ", "dsigma   ");
    table.write("# ---------", "---------");

    // compose query
    sqlitepp::statement st (session());
    st << "SELECT L,l1,l2,QUOTE(cheb) FROM 'ionf' "
            "WHERE ni = :ni "
            "  AND li = :li "
            "  AND mi = :mi "
            "  AND  S = :S  "
            "  AND Ei = :Ei "
            "ORDER BY Ei ASC",
        sqlitepp::into(L), sqlitepp::into(l1), sqlitepp::into(l2), sqlitepp::into(blob),
        sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
        sqlitepp::use(S), sqlitepp::use(Ei);

    // retrieve data (terminate if no data)
    while (st.exec())
    {
        if (not blob.empty())
        {
            // translate blob to Chebyshev coefficients
            cArray coeffs;
            coeffs.fromBlob(blob);
            Chebyshev<double,Complex> cb (coeffs, 0, 1);

            // store data
            CB[L].push_back({ l1, l2, cb });
        }
    }

    // precompute radial integrals for all combinations of angular momenta
    std::map<int,std::vector<double>> R;
    for (std::pair<int,std::vector<Info>> const & data : CB)
    {
        int L = data.first;
        std::vector<Info> const & info = data.second;

        int Nang = info.size();

        R[L].resize(Nang * Nang);

        for (int ill = 0; ill < Nang; ill++)
        for (int illp = 0; illp < Nang; illp++)
        {
            int l1 = info[ill].l1, l1p = info[illp].l1;
            int l2 = info[ill].l2, l2p = info[illp].l2;

            Complex phase0 = special::pow_int(Complex(0.,1.),-l1+l1p-l2+l2p);

            Chebyshev<double,Complex> const & f  = info[ill].cb;
            Chebyshev<double,Complex> const & fp = info[illp].cb;

            int tail = f.tail(1e-10);
            int tailp = fp.tail(1e-10);

            auto integrand = [&](double alpha) -> double
            {
                double cos_alpha = std::cos(alpha);
                double sin_alpha = std::sin(alpha);

                double k1 = kmax * cos_alpha;
                double k2 = kmax * sin_alpha;

                // relative Coulomb phase-shifts
                //   phase1 = exp(i sigma_l1(k1) - i sigma_l1p(k1))
                //   phase2 = exp(i sigma_l2(k2) - i sigma_l2p(k2))
                Complex phase1n = 1; for (int l = l1p + 1; l <= l1;  l++) phase1n *= Complex(l + 1, -1./k1);
                Complex phase1d = 1; for (int l = l1  + 1; l <= l1p; l++) phase1d *= Complex(l + 1, -1./k1);
                Complex phase2n = 1; for (int l = l2p + 1; l <= l2;  l++) phase2n *= Complex(l + 1, -1./k2);
                Complex phase2d = 1; for (int l = l2  + 1; l <= l2p; l++) phase2d *= Complex(l + 1, -1./k2);

                // asymptotic values
                if (k1 == 0) phase1n = special::pow_int(Complex(0.,-1.),l1);
                if (k1 == 0) phase1d = special::pow_int(Complex(0.,-1.),l1p);
                if (k2 == 0) phase2n = special::pow_int(Complex(0.,-1.),l2);
                if (k2 == 0) phase2d = special::pow_int(Complex(0.,-1.),l2p);

                // combined phase
                Complex phase1 = phase1n / phase1d; phase1 /= std::abs(phase1);
                Complex phase2 = phase2n / phase2d; phase2 /= std::abs(phase2);

                // radial parts (ionf) calculated by hex-ecs
                Complex u = f.clenshaw(cos_alpha, tail);
                Complex v = fp.clenshaw(cos_alpha, tailp);

                // use real part only, we are summing symmetrically, so the imaginary goes away
                Complex result = phase0 * phase1 * phase2 * u * std::conj(v);
                //std::cout << ill << " " << illp << " " << alpha << " " << result << std::endl;
                return result.real();
            };

            ClenshawCurtis<decltype(integrand),double> CC (integrand);
            R[L][ill * Nang + illp] = CC.integrate(0, special::constant::pi_quart);
        }
    }

    // process all relative angles
    for (double t : theta12)
    {
        // goniometric functions of theta
        double sint = std::sin(t * afactor);
        double cost = std::cos(t * afactor);

        // sum over all partial waves
        double dsigma = 0;

        // process all partial wave contributions
        for (std::pair<int,std::vector<Info>> const & data : CB)
        {
            int L = data.first;
            std::vector<Info> const & info = data.second;

            int Nang = info.size();

            for (int ill = 0; ill < Nang; ill++)
            for (int illp = 0; illp < Nang; illp++)
            {
                int l1 = info[ill].l1, l1p = info[illp].l1;
                int l2 = info[ill].l2, l2p = info[illp].l2;

                int minlambda = std::max(std::abs(l1-l1p),std::abs(l2-l2p));
                int maxlambda = std::min(         l1+l1p ,         l2+l2p );

                for (int lambda = minlambda; lambda <= maxlambda; lambda++)
                {
                    double angfactor = (2*lambda+1) * special::pow_int(-1, lambda)
                        * special::computef(lambda,l1,l2,l1p,l2p,L)
                        * 0.5 * gsl_sf_legendre_Pl(lambda, cost); // * sint;

                    //std::cout << "t = " << t << ", L = " << L << ", ill = " << ill << ", illp = " << illp << ", lambda = " << lambda << ", angfactor = " << angfactor << ", R = " << R[L][ill*Nang+illp] << std::endl;
                    dsigma += R[L][ill * Nang + illp] * angfactor;
                }
            }
        }

        // write line to table
        table.write(t, dsigma * lfactor * lfactor * (2*S + 1) / 4 / std::sqrt(Ei));
    }

    return true;
}
