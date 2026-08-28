#ifdef CG_BLOCK_JACOBI

    // diagonal block matrices
    std::vector<CsrMatrix<LU_int_t,Complex>> blockmat (N);
    
    // diagonal block LU decompositions
    std::vector<std::shared_ptr<LUft>> blocklu (N);

    // prepare all diagonal block decompositions
    for (int i = 0; i < N; i++)
    {
        blockmat[i] = (bD(i,i) * bS + bS(i,i) * bD).tocoo<LU_int_t>().tocsr();
        
        blocklu[i].reset(LUft::Choose("umfpack"));
        blocklu[i]->factorize(blockmat[i]);
    }

#endif
