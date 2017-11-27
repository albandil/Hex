
    // evaluate source/solution fields in selected points
    rArray pts = linspace(bspline.Rmin(), bspline.Rmax(), 1001);
    cArray evalS = Bspline::zip(bspline, bspline, S.data(), pts, pts);
    cArray evalT = Bspline::zip(bspline, bspline, T.data(), pts, pts);
    
    // write source/solution as VTK datasets
    VTKRectGridFile vtk;
    vtk.setGridX(pts);
    vtk.setGridY(pts);
    vtk.setGridZ(rArray{0});
    vtk.appendScalarAttribute("S", realpart(evalS));
    vtk.appendScalarAttribute("T", realpart(evalT));
    vtk.writePoints("fields.vtk");

    return 0;
