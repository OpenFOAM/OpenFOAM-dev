{
    const volScalarField& psi = thermo.psi();

    tmp<volVectorField> tHbyA;
    if (pressureImplicitPorosity)
    {
        tHbyA = constrainHbyA(trTU()&UEqn.H(), U, p);
    }
    else
    {
        tHbyA = constrainHbyA(trAU()*UEqn.H(), U, p);
    }
    volVectorField& HbyA = tHbyA.ref();

    tUEqn.clear();

    bool closedVolume = false;

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    );
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    closedVolume = adjustPhi(phiHbyA, U, p);

    while (simple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;

        if (pressureImplicitPorosity)
        {
            tpEqn =
            (
                fvm::laplacian(rho*trTU(), p)
              + fvModels.source(psi, p, rho.name())
             ==
                fvc::div(phiHbyA)
            );
        }
        else
        {
            tpEqn =
            (
                fvm::laplacian(rho*trAU(), p)
              + fvModels.source(psi, p, rho.name())
             ==
                fvc::div(phiHbyA)
            );
        }

        fvScalarMatrix& pEqn = tpEqn.ref();

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    if (pressureImplicitPorosity)
    {
        U = HbyA - (trTU() & fvc::grad(p));
    }
    else
    {
        U = HbyA - trAU()*fvc::grad(p);
    }

    U.correctBoundaryConditions();
    fvConstraints.constrain(U);

    pressureControl.limit(p);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();
    }

    rho = thermo.rho();
    rho.relax();
}
