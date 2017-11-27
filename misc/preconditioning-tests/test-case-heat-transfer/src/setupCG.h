
    std::array<std::string,9> gauge =
    {{
        "\x1b[48;5;21m\x1b[38;5;196m \e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▁\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▂\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▃\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▄\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▅\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▆\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m▇\e[0m",
        "\x1b[48;5;21m\x1b[38;5;196m█\e[0m"
    }};
    
    RowMatrix<Complex> E (N, N), F (N, N), G (N, N);
    
    // sort 2D eigenvalues
    std::vector<std::pair<int,int>> e2d (N*N);
    for (unsigned i = 0; i < e2d.size(); i++)
    {
        e2d[i] = std::make_pair(i / N, i % N);
    }
    std::sort
    (
        e2d.begin(), e2d.end(),
        [&](std::pair<int,int> const & u, std::pair<int,int> const & v) -> bool
        {
            return e[u.first].real() + e[u.second].real()
                 < e[v.first].real() + e[v.second].real();
        }
    );
    
    ConjugateGradients<Complex,cArray,cArrayView> CG;
    CG.reset();
    CG.verbose              = true;
    CG.stationary           = stationary;
    CG.new_array            = [&](std::size_t n, std::string s) -> cArray { return cArray(n); };
    CG.compute_norm         = [&](cArrayView r) -> Real  { return r.norm(); };
    CG.scalar_product       = [&](const cArrayView x, const cArrayView y) -> Complex { return x | y; };
    CG.axby                 = [&](Complex a, cArrayView x, Complex b, cArrayView y) { x = a * x + b * y; };
    CG.matrix_multiply      = [&](const cArrayView p, cArrayView q) { q = A.dot(p); };
    CG.process_solution     = [&](int iter, const cArrayView x)
    {
        RowMatrix<Complex> R (N, N, S.data() - A.dot(x));
        
        // extract eigen-components of the right-hand side
        blas::gemm(1.0, S,     C, 0.0, E);
        blas::gemm(1.0, C.T(), E, 0.0, F);
        
        // extract eigen-components of the residual
        blas::gemm(1.0, R,     C, 0.0, E);
        blas::gemm(1.0, C.T(), E, 0.0, G);
        
        // resample
        int m = N;
        rArray eb (m), er (m);
        for (int i = 0; i < N*N; i++)
        {
            eb[i * m / (N*N)] += sqrabs(F(e2d[i].first, e2d[i].second));
            er[i * m / (N*N)] += sqrabs(G(e2d[i].first, e2d[i].second));
        }
        
        // write a graphical representation
        std::cout << "| ";
        for (int k = 0; k < m; k++)
        {
            int val = std::max(0., std::min(8., std::log10(1e8 * std::sqrt(er[k] / eb[k]))));
            std::cout << gauge[val];
        }
    };
    CG.apply_preconditioner = [&](const cArrayView r, cArrayView z)
