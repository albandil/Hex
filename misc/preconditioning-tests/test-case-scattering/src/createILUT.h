#ifdef CG_ILUT
    // calculate incomplete LU decomposition
    std::shared_ptr<LUft> ilut (LUft::Choose("umfpack"));
    ilut->rdata("drop_tolerance") = droptol;
    ilut->factorize(A);
    ilut->get().plot(format("lu-%f.png", droptol));
    std::cout << "ILUT(" << droptol << ") nonzeros: " << ilut->size()/16 - N*N << "\n\n";
#endif
