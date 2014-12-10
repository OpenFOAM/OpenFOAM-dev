{
    volScalarField thermoRho = psi*p + alphal*rhol0;

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
    (
        fvc::domainIntegrate(mag(rho - thermoRho))/totalMass
    ).value();

    scalar globalContErr =
    (
        fvc::domainIntegrate(rho - thermoRho)/totalMass
    ).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
}
