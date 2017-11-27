{
    // 2D eigenvalues
    ColMatrix<Complex> ca = A.torow().T();
    ColMatrix<Complex> cs = kron(bS,bS).tocoo<LU_int_t>().torow().T();
    cArray ee (N * N);
    ca.diagonalize_g(cs, ee, nullptr, nullptr);
    rArray rr = realpart(ee);
    std::sort(rr.begin(), rr.end());
    write_array(rr, "2D-eigenvalues.txt");
}
