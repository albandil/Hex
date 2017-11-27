
    // solution and right-hand side (sources)
    RowMatrix<Complex> S (N,N), T (N,N), T0 (N,N);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        Real r1 = 0.5 * (bspline.t(i).real() + bspline.t(i + 1).real());
        Real r2 = 0.5 * (bspline.t(j).real() + bspline.t(j + 1).real());
        S(i,j) = r1 * r2 * std::exp(-r1) * std::exp(-r2);
    }
