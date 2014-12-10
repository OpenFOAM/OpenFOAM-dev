{
    volScalarField rAUrel(1.0/UrelEqn().A());
    volVectorField HbyA("HbyA", Urel);
    HbyA = rAUrel*UrelEqn().H();
    UrelEqn.clear();

    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf());
    adjustPhi(phiHbyA, Urel, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUrel, p) == fvc::div(phiHbyA)
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
    Urel = HbyA - rAUrel*fvc::grad(p);
    Urel.correctBoundaryConditions();
    fvOptions.correct(Urel);
}
