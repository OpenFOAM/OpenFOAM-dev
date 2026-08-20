{
    volInternalScalarField conserve(-fvi::div(mesh.phi()));
    // The ddt term constructed by hand because it would be wrong for
    // Backward Differencing in time.

    conserve += (1.0 - mesh.V0()/mesh.V())/runTime.deltaTValue();

    scalar sumLocalContErr = runTime.deltaTValue()*
        weightedAverage(eval(mag(conserve))(), mesh.V()).value();

    scalar globalContErr = runTime.deltaTValue()*
        weightedAverage(conserve, mesh.V()).value();

    Info<< "volume continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr << endl;
}
