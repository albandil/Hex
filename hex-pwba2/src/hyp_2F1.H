#ifndef HYP_2F1_H
#define HYP_2F1_H

#include <complex>

std::complex<double> hyp_2F1
(
    std::complex<double> const & a,
    std::complex<double> const & b,
    std::complex<double> const & c,
    std::complex<double> const & z
);

#endif