    std::cout << "\n   iteration | time        | residual       | eigensubspace populations" << std::endl;
    CG.solve(chi, psi, 1e-5, 0, 1000);
