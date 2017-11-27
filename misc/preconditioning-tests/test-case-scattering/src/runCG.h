    std::cout << "\n   iteration | time        | residual       | eigensubspace populations" << std::endl;
    CG.solve(chi.data(), psi.data(), 1e-5, 0, 1000);
