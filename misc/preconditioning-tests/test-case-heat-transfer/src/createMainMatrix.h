
    // matrix of the 2D system
    CsrMatrix<LU_int_t,Complex> A = (kron(bD,bS) + kron(bS,bD)).tocoo<LU_int_t>().tocsr();
    A.plot("A.png");
    std::cout << "\nA   nonzeros: " << A.p().back() << std::endl;
