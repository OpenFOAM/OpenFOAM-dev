Info<< "magnetic flux divergence error = "
    << runTime.deltaTValue()
      *mag(fvc::div(phiB))().weightedAverage(mesh.V()).value()
    << endl;
