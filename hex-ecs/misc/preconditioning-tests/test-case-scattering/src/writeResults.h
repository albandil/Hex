    
    // evaluate source/solution fields in selected points
    rArray pts = linspace(bspline.Rmin(), bspline.Rmax(), 1001);
    cArray evalChi = Bspline::zip(bspline, bspline, chi.data(), pts, pts);
    cArray evalPsi = Bspline::zip(bspline, bspline, psi.data(), pts, pts);
    
    // write source/solution as VTK datasets
    VTKRectGridFile vtk;
    vtk.setGridX(pts);
    vtk.setGridY(pts);
    vtk.setGridZ(rArray{0});
    vtk.appendScalarAttribute("Re-chi", realpart(evalChi));
    vtk.appendScalarAttribute("Im-chi", imagpart(evalChi));
    vtk.appendScalarAttribute("Re-psi", realpart(evalPsi));
    vtk.appendScalarAttribute("Im-psi", imagpart(evalPsi));
    vtk.writePoints("fields.vtk");
    
    return 0;
