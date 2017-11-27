
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
    
    std::cout << "\n   iteration | time        | residual       | eigensubspace populations" << std::endl;
    
    Timer timer;
    
    RowMatrix<Complex> R (N, N, S.data() - A.dot(T));
    
    for (int k = 1; k < 1000; k++)
    {
        std::cout << '\t';
        std::cout << std::setw(4) << std::right << k;
        std::cout << " | ";
        std::cout << std::setw(11) << std::left << timer.nice_time();
        std::cout << " | ";
        std::cout << std::setw(15) << std::left << R.data().norm() / S.data().norm();
        
        // stationary iteration
        
        
        // update residual
        R = RowMatrix<Complex>(N, N, S.data() - A.dot(x));
        
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
    }
