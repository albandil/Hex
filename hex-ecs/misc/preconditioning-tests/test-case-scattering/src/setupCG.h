
    ConjugateGradients<Complex,cArray,cArrayView> CG;
    CG.reset();
    CG.verbose              = true;
    CG.stationary           = stationary;
    CG.new_array            = [&](std::size_t n, std::string s) -> cArray { return cArray(n); };
    CG.compute_norm         = [&](cArrayView r) -> double  { return r.norm(); };
    CG.scalar_product       = [&](const cArrayView x, const cArrayView y) -> Complex { return x | y; };
    CG.axby                 = [&](Complex a, cArrayView x, Complex b, cArrayView y) { x = a * x + b * y; };
    CG.matrix_multiply      = [&](const cArrayView p, cArrayView q) { q = A.dot(p); };
    CG.process_solution     =  [&](unsigned iter, const cArrayView x)
    {
        std::cout << "| ";
    };
    CG.apply_preconditioner = [&](const cArrayView r, cArrayView z)
