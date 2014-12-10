{
    volVectorField HbyA("HbyA", Uc);
    HbyA = rAUc*UcEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           (fvc::interpolate(HbyA) & mesh.Sf())
         + alphacf*rAUcf*fvc::ddtCorr(Uc, phic)
         + phicForces
        )
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p.boundaryField(),
        (
            phiHbyA.boundaryField()
            - (mesh.Sf().boundaryField() & Uc.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUcf.boundaryField())
    );

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rAUcf, p)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - pEqn.flux()/alphacf;

            p.relax();

            Uc = HbyA
              + rAUc*fvc::reconstruct((phicForces - pEqn.flux()/alphacf)/rAUcf);
            Uc.correctBoundaryConditions();
        }
    }
}

#include "continuityErrs.H"
