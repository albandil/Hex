#ifdef CG_DILU

    CsrMatrix<LU_int_t,Complex> dilu (A);
    
    // pivots
    cArray d = A.diag();
    
    // loop over all matrix rows
    for (int i = 0; i < N*N; i++)
    {
        // loop over all non-zero elements in a row
        for (int x = dilu.p()[i]; x < dilu.p()[i + 1]; x++)
        {
            // current column
            int j = dilu.i()[x];
            
            // only do something on diagonal
            if (i == j)
            {
                // update current diagonal element
                dilu.x()[x] = 1.0 / d[i];
                
                // update all further pivots
                for (int y = x + 1; y < dilu.p()[i + 1]; y++)
                {
                    // other column
                    int k = dilu.i()[y];
                    
                    // off-diagonal element (aik = aki)
                    Complex aik = dilu.x()[y];
                    
                    // update graph-connected pivot
                    d[j] -= aik * aik / d[i];
                }
            }
        }
    }
    
    write_array(realpart(d), "dilu-inverse-pivots.txt");

#endif
