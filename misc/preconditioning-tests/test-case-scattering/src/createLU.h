
    // calculate exact LU decomposition
    std::shared_ptr<LUft> lu0;
    lu0.reset(LUft::Choose("umfpack"));
    lu0->factorize(A);
    lu0->solve(chi.data(), psi0.data(), 1);
    std::cout << "LU  nonzeros: " << lu0->size()/16 - N*N << std::endl;
