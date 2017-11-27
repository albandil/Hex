{
        ColMatrix<Complex> ca = A.torow().T();
        ColMatrix<Complex> cs = kron(bS,bS).tocoo<LU_int_t>().torow().T();
        cArray ee (N * N);
        ca.diagonalize_g(cs, ee, nullptr, nullptr);
        cArray rr = -ee;
        std::sort(rr.begin(), rr.end(), Complex_realpart_less);
        write_array(rr, "2D-eigenvalues.txt");
}
