{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA)
     ==
        parcels.Srho()
      + fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            pDDtEqn
          - fvm::laplacian(rhorAUf, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }

    p.relax();

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H" // NOTE: flux and time scales now inconsistent
    #include "compressibleContinuityErrs.H"

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);

    Info<< "p min/max = " << min(p).value() << ", " << max(p).value() << endl;
}
