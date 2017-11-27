
    // B-spline basis
    int order = 4;
    Real theta = 0;
    Bspline bspline (order, theta, rArray{}, knots, rArray{ knots.back() });
    int N = bspline.Nspline();
    std::cout << "B-spline knots: " << knots << std::endl;
    std::cout << "\nB-spline count: " << N << std::endl;
