
    // matrix of the 3D system
    CsrMatrix<LU_int_t,Complex> A;
    {
        std::size_t vol = (2 * order + 1) * N;
        std::size_t dim = N*N*N;
        
        NumberArray<LU_int_t> I (vol*vol*vol);
        NumberArray<LU_int_t> J (vol*vol*vol);
        NumberArray<Complex>  V (vol*vol*vol);
        
        # pragma omp parallel for
        
        for (int xi = 0; xi < N; xi++)
        for (int xj = 0; xj < N; xj++)
        if (std::max(xi,xj) <= std::min(xi,xj) + order)
        
        for (int yi = 0; yi < N; yi++)
        for (int yj = 0; yj < N; yj++)
        if (std::max(yi,yj) <= std::min(yi,yj) + order)
        
        for (int zi = 0; zi < N; zi++)
        for (int zj = 0; zj < N; zj++)
        if (std::max(zi,zj) <= std::min(zi,zj) + order)
        
        {
            std::size_t X = xi * (2 * order + 1) + order + xi - xj;
            std::size_t Y = yi * (2 * order + 1) + order + yi - yj;
            std::size_t Z = zi * (2 * order + 1) + order + zi - zj;
            
            I[(X * vol + Y) * vol + Z] = (xi * N + yi) * N + zi;
            J[(X * vol + Y) * vol + Z] = (xj * N + yj) * N + zj;
            V[(X * vol + Y) * vol + Z] =
            (
                Etot * bS(xi,xj) * bS(yi,yj) * bS(zi,zj)
                - (0.5 * bD(xi,xj) - 2.0 * bM(xi,xj)) * bS(yi,yj) * bS(zi,zj)
                - (0.5 * bD(yi,yj) - 2.0 * bM(yi,yj)) * bS(zi,zj) * bS(xi,xj)
                - (0.5 * bD(zi,zj) - 2.0 * bM(zi,zj)) * bS(xi,xj) * bS(yi,yj)
                - bR(xi,yi,xj,yj) * bS(zi,zj)
                - bR(yi,zi,yj,zj) * bS(xi,xj)
                - bR(zi,xi,zj,xj) * bS(yi,yj)
            );
        }
        
        A = CooMatrix<LU_int_t,Complex>(dim, dim, std::move(I), std::move(J), std::move(V)).tocsr();
    }
    //A.plot("A.png");
    std::cout << "\nA   nonzeros: " << A.p().back() << std::endl;
    
