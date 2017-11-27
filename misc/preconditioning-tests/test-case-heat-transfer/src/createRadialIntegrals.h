
    // B-spline integrals
    RadialIntegrals rad (bspline, bspline, 0);
    rad.verbose(false);
    rad.setupOneElectronIntegrals();
    SymBandMatrix<Complex> const & bS = rad.S();
    SymBandMatrix<Complex> const & bD = rad.D();
    bS.tocoo<LU_int_t>().tocsr().plot("S.png");
    bD.tocoo<LU_int_t>().tocsr().plot("D.png");
