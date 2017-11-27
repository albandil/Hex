    // B-spline integrals
    RadialIntegrals rad (bspline, bspline, 1);
    rad.verbose(false);
    rad.setupOneElectronIntegrals(par, cmd);
    rad.setupTwoElectronIntegrals(par, cmd);
    SymBandMatrix<Complex> const & bS = rad.S();
    SymBandMatrix<Complex> const & bD = rad.D();
    SymBandMatrix<Complex> const & bM = rad.Mm1();
    BlockSymBandMatrix<Complex> const & bR = rad.R_tr_dia(0);
    bS.tocoo<LU_int_t>().tocsr().plot("S.png");
    bD.tocoo<LU_int_t>().tocsr().plot("D.png");
