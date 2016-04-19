{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    tUEqn.clear();

    bool closedVolume = false;

    if (simple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi)
           *fvc::flux(HbyA)
        );

        MRF.makeRelative(fvc::interpolate(psi), phid);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::div(phid, p)
              - fvm::laplacian(rhorAUf, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi == pEqn.flux();
            }
        }
    }
    else
    {
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(rho*HbyA));
        MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

        closedVolume = adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA)
              - fvm::laplacian(rhorAUf, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }


    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);

    if (!simple.transonic())
    {
        rho.relax();
    }

    Info<< "rho max/min : "
        << max(rho).value() << " "
        << min(rho).value() << endl;
}
