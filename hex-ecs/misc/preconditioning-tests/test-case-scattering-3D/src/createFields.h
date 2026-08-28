
    // solution and right-hand side (sources)
    cArray chi (N*N*N), psi (N*N*N), psi0 (N*N*N);
    for (int i = 0; i < N - 1; i++)
    for (int j = 0; j < N - 1; j++)
    for (int k = 0; k < N - 1; k++)
    {
        Real r1 = 0.5 * (bspline.t(i).real() + bspline.t(i + 1).real());
        Real r2 = 0.5 * (bspline.t(j).real() + bspline.t(j + 1).real());
        Real r3 = 0.5 * (bspline.t(k).real() + bspline.t(k + 1).real());
        chi[(i * N + j) * N + k] = r1 * r2 * r3 * std::exp(-r1) * std::exp(-r2) * std::exp(-r3);
    }
