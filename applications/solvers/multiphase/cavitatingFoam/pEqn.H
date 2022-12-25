{
    if (pimple.nCorr() == 1)
    {
        mixture.correctPressure();
    }

    surfaceScalarField rhof("rhof", fvc::interpolate(rho));

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    phi = fvc::flux(HbyA) + rhorAUf*fvc::ddtCorr(U, phi, Uf);
    fvc::makeRelative(phi, U);

    surfaceScalarField phiGradp(rhorAUf*mesh.magSf()*fvc::snGrad(p));

    phi -= phiGradp/rhof;

    const volScalarField rho0(rho - mixture.psi()*p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho)
          + mixture.psi()*correction(fvm::ddt(p))
          + fvc::div(phi, rho)
          + fvc::div(phiGradp)
          - fvm::laplacian(rhorAUf, p)
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi += (phiGradp + pEqn.flux())/rhof;
        }
    }

    Info<< "Predicted p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    rho == rho0 + mixture.psi()*p;

    mixture.correct();

    mixture.correctPressure();

    Info<< "Phase-change corrected p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    // Correct velocity

    U = HbyA - rAU*fvc::grad(p);

    // Remove the swirl component of velocity for "wedge" cases
    if (pimple.dict().found("removeSwirl"))
    {
        label swirlCmpt(pimple.dict().lookup<label>("removeSwirl"));

        Info<< "Removing swirl component-" << swirlCmpt << " of U" << endl;
        U.field().replace(swirlCmpt, 0.0);
    }

    U.correctBoundaryConditions();

    Info<< "max(U) " << max(mag(U)).value() << endl;

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, fvc::absolute(phi, U));
}
