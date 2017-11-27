    // B-spline basis
    int order = 4;
    Real theta = 0.63;
    Bspline bspline (order, theta, rArray{}, rknots, cknots);
    int N = bspline.Nspline();
    std::cout << "B-spline rknots: " << bspline.rknots() << std::endl;
    std::cout << "\nB-spline cknots: " << bspline.cknots2() << std::endl;
    std::cout << "\nB-spline count: " << N << std::endl;
    
