{
    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn().A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    UEqn.clear();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*(fvc::interpolate(HbyA) & mesh.Sf())
    );

    fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rho*rAU, p)
         ==
            parcels.Srho()
          + fvOptions(psi, p, rho.name())
        );

        fvOptions.constrain(pEqn);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }

    p.relax();

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "compressibleContinuityErrs.H"

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);
    rho.relax();

    Info<< "p min/max = " << min(p).value() << ", " << max(p).value() << endl;
}
