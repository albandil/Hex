
    // matrix of the 2D system
    CsrMatrix<LU_int_t,Complex> A;
    {
        std::size_t vol = (2 * order + 1) * N;
        
        NumberArray<LU_int_t> I;  I.reserve(vol*vol);
        NumberArray<LU_int_t> J;  J.reserve(vol*vol);
        NumberArray<Complex>  V;  V.reserve(vol*vol);
        
        for (int xi = 0; xi < N; xi++)
        for (int xj = 0; xj < N; xj++)
        if (std::max(xi,xj) <= std::min(xi,xj) + order)
        
        for (int yi = 0; yi < N; yi++)
        for (int yj = 0; yj < N; yj++)
        if (std::max(yi,yj) <= std::min(yi,yj) + order)
        
        {
            I.push_back(xi * N + yi);
            J.push_back(xj * N + yj);
            V.push_back
            (
                Etot * bS(xi,xj) * bS(yi,yj)
                - (0.5 * bD(xi,xj) - bM(xi,xj)) * bS(yi,yj)
                - (0.5 * bD(yi,yj) - bM(yi,yj)) * bS(xi,xj)
                - bR(xi,yi,xj,yj)
            );
        }
        
        A = CooMatrix<LU_int_t,Complex>(N*N, N*N, I, J, V).tocsr();
    }
    A.plot("A.png");
    std::cout << "\nA   nonzeros: " << A.p().back() << std::endl;
    
