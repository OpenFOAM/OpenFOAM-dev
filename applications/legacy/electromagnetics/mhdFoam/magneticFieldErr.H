Info<< "magnetic flux divergence error = "
    << runTime.deltaTValue()
      *weightedAverage(eval(mag(fvi::div(phiB)))(), mesh.V()).value()
    << endl;
