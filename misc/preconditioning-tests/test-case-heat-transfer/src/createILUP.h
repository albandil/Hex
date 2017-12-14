#ifdef CG_ILUP
    // first construct adjacency graph of the matrix
    std::vector<std::vector<int>> adjc (N*N), adjl (N*N);
    for (int row = 0; row < N*N; row++)
    for (int idx = A.p()[row]; idx < A.p()[row + 1]; idx++)
    {
        int col = A.i()[idx];
        if (row <= col)
        {
            adjc[row].push_back(col);
            adjl[row].push_back(0);
        }
    }
    
    // expand the adjacency graph to the given level
    for (int row = 0; row < N*N; row++)
    {
        # pragma omp parallel for
        for (unsigned i = 0; i < adjc[row].size(); i++)
        {
            int nbr = adjc[row][i];
            int lvl = adjl[row][i];
            
            for (unsigned j = 0; j < adjc[row].size(); j++)
            {
                int conn = adjc[row][j];
                int clvl = adjl[row][j];
                
                if (conn > nbr and std::max(lvl, clvl) < level)
                {
                    std::size_t pos = std::lower_bound(adjc[nbr].begin(), adjc[nbr].end(), conn) - adjc[nbr].begin();
                    if (pos < adjc[nbr].size() and adjc[nbr][pos] == conn)
                    {
                        adjl[nbr][pos] = std::min(adjl[nbr][pos], std::max(lvl, clvl) + 1);
                    }
                    if (pos == adjc[nbr].size() or adjc[nbr][pos] != conn)
                    {
                        adjc[nbr].insert(adjc[nbr].begin() + pos, conn);
                        adjl[nbr].insert(adjl[nbr].begin() + pos, std::max(lvl, clvl) + 1);
                    }
                }
            }
        }
    }
    
    // print statistics
    std::vector<int> stats = { N*N };
    for (int row = 0; row < N*N; row++)
    for (unsigned i = 0; i < adjc[row].size(); i++)
    {
        int col = adjc[row][i];
        int lvl = adjl[row][i];
        
        stats.resize(std::max<int>(stats.size(), lvl + 1));
        
        if (row != col)
            stats[lvl] += 2;
    }
    std::size_t nnz = std::accumulate(stats.begin(), stats.end(), 0);
    std::cout << "\nILUP(" << level << ") nonzeros: " << nnz << std::endl;
    for (unsigned lvl = 0; lvl < stats.size(); lvl++)
    {
        std::cout << "  - level " << lvl << ": " << stats[lvl] << std::endl;
    }
    std::vector<std::vector<int>>().swap(adjl);
    
    // allocate CSR matrix for the L' factor
    NumberArray<LU_int_t> P (N*N + 1);
    NumberArray<LU_int_t> I ((nnz + N*N)/2);
    NumberArray<Complex>  X ((nnz + N*N)/2);
    for (int i = 0, n = 0; i < N*N; i++)
    {
        int idx = A.p()[i];
        
        for (int j : adjc[i])
        {
            I[n] = j;
            
            while (idx < A.p()[i + 1] and A.i()[idx] < j)
                idx++;
            
            if (A.i()[idx] == j)
                X[n] = A.x()[idx];
            
            n++;
        }
        
        P[i + 1] = n;
    }
    
    // diagonal of the LDL' decomposition
    cArray ilupd (N*N);
    
    // calculate the L' factor
    for (int k = 0; k < N*N; k++)
    {
        // backup pivot
        ilupd[k] = X[P[k]];
        
        // loop over all rows that have non-zero in k-th column (skip self)
        # pragma omp parallel for
        for (LU_int_t idx1 = P[k]; idx1 < P[k + 1]; idx1++) if (I[idx1] != k)
        {
            // row index
            int i = I[idx1];
            
            // loop over all entries (columns) of the k-th and i-th row, synchronously
            LU_int_t idx2 = P[k], idx3 = P[i];
            while (idx2 < P[k + 1] and idx3 < P[i + 1])
            {
                // column index of entry on k-th row
                int j = I[idx2];
                
                // skip pivot (it just trivially annihilates the column)
                if (j != k)
                {
                    // check if the i-th row has also the column j by advancing position to the first non-zero entry with not-less column index
                    while (idx3 < P[i + 1] and I[idx3] < j)
                        idx3++;
                    
                    // update the element or pivot (if using MILU(p))
                    if (I[idx3] == j)
                        X[idx3] -= X[idx1] * X[idx2] / ilupd[k];
                    else if (modified)
                        X[P[i]] += std::abs(X[idx1] * X[idx2] / ilupd[k]);
                }
                
                // move to the next element of the k-th row
                idx2++;
            }
        }
        
        // ... and finally calculate L column
        for (LU_int_t idx = P[k]; idx < P[k + 1]; idx++)
        {
            X[idx] /= ilupd[k];
        }
    }
    
    // create the L matrix
    CsrMatrix<LU_int_t,Complex> ilupU (std::size_t(N*N), std::size_t(N*N), std::move(P), std::move(I), std::move(X));
    CsrMatrix<LU_int_t,Complex> ilupL = ilupU.tocoo().transpose().tocsr();
    write_array(ilupd, "ILUP-" + std::to_string(level) + ".txt");
#endif
