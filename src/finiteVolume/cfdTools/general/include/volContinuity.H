{
    volScalarField conserve(-fvc::div(mesh.phi()));
    // The ddt term constructed by hand because it would be wrong for
    // Backward Differencing in time.

    conserve.primitiveFieldRef() +=
        (1.0 - mesh.V0()/mesh.V())/runTime.deltaTValue();

    scalar sumLocalContErr = runTime.deltaTValue()*
        mag(conserve)().weightedAverage(mesh.V()).value();

    scalar globalContErr = runTime.deltaTValue()*
        conserve.weightedAverage(mesh.V()).value();

    Info<< "volume continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr << endl;
}
