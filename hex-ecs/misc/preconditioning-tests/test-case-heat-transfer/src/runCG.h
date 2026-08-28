
    std::cout << "\n   iteration | time        | residual       | eigensubspace populations" << std::endl;
    CG.solve(S.data(), T.data(), 1e-5, 0, 1000);
