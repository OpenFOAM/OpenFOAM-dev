{
    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));
    volVectorField HbyASp(rASpUc/rAUc*HbyA);

    surfaceScalarField phiHbyASp
    (
        "phiHbyASp",
        (
           fvc::flux(HbyASp)
         + alphacf*rASpUcf*fvc::ddtCorr(Uc, phic, Ucf)
        )
    );

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyASp, Uc);
        adjustPhi(phiHbyASp, Uc, p);
        fvc::makeAbsolute(phiHbyASp, Uc);
    }

    phiHbyASp += phicSUSu;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Uc, phiHbyASp, rASpUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rASpUcf, p)
         ==
            fvc::ddt(alphac)
          + fvc::div(alphacf*phiHbyASp)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyASp - pEqn.flux()/alphacf;

            // Explicitly relax pressure for momentum corrector
            p.relax();

            Uc =
                HbyA
              + rAUc
               *fvc::reconstruct
                (
                    (phicSUSu + phicSUSp*phic - pEqn.flux()/alphacf)/rASpUcf
                );
            Uc.correctBoundaryConditions();
            fvOptions.correct(Uc);

            // Correct Ucf if the mesh is moving
            fvc::correctUf(Ucf, Uc, phic);

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phic, Uc);
        }
    }
}

#include "continuityErrs.H"
