    // calculate exact LU decomposition
    std::shared_ptr<LUft> lu0;
    std::cout << "Direct solution (exact)\n" << std::endl;
    lu0.reset(LUft::Choose("umfpack"));
    lu0->factorize(A);
    lu0->solve(chi, psi0, 1);
