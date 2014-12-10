{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
    (
        fvc::domainIntegrate(mag(rho - thermo.rho()))/totalMass
    ).value();

    scalar globalContErr =
    (
        fvc::domainIntegrate(rho - thermo.rho())/totalMass
    ).value();

    cumulativeContErr[i] += globalContErr;

    Info<< "time step continuity errors (" << mesh.name() << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr[i]
        << endl;
}
