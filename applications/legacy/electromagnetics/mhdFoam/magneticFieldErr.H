Info<< "magnetic flux divergence error = "
    << runTime.deltaTValue()
      *weightedAverage(mag(fvi::div(phiB)), mesh.V()).value()
    << endl;
