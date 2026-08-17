Info<< "magnetic flux divergence error = "
    << runTime.deltaTValue()
      *mag(fvi::div(phiB))().weightedAverage(mesh.V()).value()
    << endl;
