{
    if (pimple.nCorrPIMPLE() == 1)
    {
        p =
        (
            rho
          - alphal*rhol0
          - ((alphav*psiv + alphal*psil) - psi)*pSat
        )/psi;
    }

    surfaceScalarField rhof("rhof", fvc::interpolate(rho));

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    phi = (fvc::interpolate(HbyA) & mesh.Sf())
         + rhorAUf*fvc::ddtCorr(U, phi);

    surfaceScalarField phiGradp(rhorAUf*mesh.magSf()*fvc::snGrad(p));

    phi -= phiGradp/rhof;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          - (rhol0 + (psil - psiv)*pSat)*fvc::ddt(alphav) - pSat*fvc::ddt(psi)
          + fvc::div(phi, rho)
          + fvc::div(phiGradp)
          - fvm::laplacian(rhorAUf, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += (phiGradp + pEqn.flux())/rhof;
        }
    }

    Info<< "Predicted p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    rho == max
    (
        psi*p
      + alphal*rhol0
      + ((alphav*psiv + alphal*psil) - psi)*pSat,
        rhoMin
    );

    #include "alphavPsi.H"

    p =
    (
        rho
      - alphal*rhol0
      - ((alphav*psiv + alphal*psil) - psi)*pSat
    )/psi;

    p.correctBoundaryConditions();

    Info<< "Phase-change corrected p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    // Correct velocity

    U = HbyA - rAU*fvc::grad(p);

    // Remove the swirl component of velocity for "wedge" cases
    if (pimple.dict().found("removeSwirl"))
    {
        label swirlCmpt(readLabel(pimple.dict().lookup("removeSwirl")));

        Info<< "Removing swirl component-" << swirlCmpt << " of U" << endl;
        U.field().replace(swirlCmpt, 0.0);
    }

    U.correctBoundaryConditions();

    Info<< "max(U) " << max(mag(U)).value() << endl;
}
