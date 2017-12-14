#ifdef CG_BLOCK_JACOBI

    // diagonal block matrices
    std::vector<CsrMatrix<LU_int_t,Complex>> blockmat (N);
    
    // diagonal block LU decompositions
    std::vector<std::shared_ptr<LUft>> blocklu (N);

    // prepare all diagonal block decompositions
    for (int i = 0; i < N; i++)
    {
        blockmat[i] = (
            Etot * bS(i,i) * bS
                - (0.5 * bD(i,i) - bM(i,i)) * bS
                - bS(i,i) * (0.5_z * bD - bM)
                - bR(i,i)
        ).tocoo<LU_int_t>().tocsr();
        
        blocklu[i].reset(LUft::Choose("umfpack"));
        blocklu[i]->factorize(blockmat[i]);
    }

#endif
