{
    rAU = 1.0/UEqn().A();
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    if (pimple.nCorrPISO() <= 1)
    {
        UEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + rAUf*fvc::ddtCorr(U, Uf)
    );

    if (p_gh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_gh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_gh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_ghEqn
        (
            fvm::laplacian(rAUf, p_gh) == fvc::div(phiHbyA)
        );

        p_ghEqn.setReference(p_ghRefCell, p_ghRefValue);

        p_ghEqn.solve(mesh.solver(p_gh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_ghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_gh.relax();

            U = HbyA - rAU*fvc::grad(p_gh);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p = p_gh + (g & mesh.C());
}
