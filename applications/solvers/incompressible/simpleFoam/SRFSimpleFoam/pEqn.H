{
    volScalarField rAUrel(1.0/UrelEqn.A());
    volVectorField HbyA("HbyA", Urel);
    HbyA = rAUrel*UrelEqn.H();

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, Urel, p);

    tmp<volScalarField> rAtUrel(rAUrel);

    if (simple.consistent())
    {
        rAtUrel = 1.0/(1.0/rAUrel - UrelEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtUrel() - rAUrel)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAUrel - rAtUrel())*fvc::grad(p);
    }

    tUrelEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Urel, phiHbyA, rAtUrel());

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtUrel(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    Urel = HbyA - rAtUrel()*fvc::grad(p);
    Urel.correctBoundaryConditions();
    fvOptions.correct(Urel);
}
