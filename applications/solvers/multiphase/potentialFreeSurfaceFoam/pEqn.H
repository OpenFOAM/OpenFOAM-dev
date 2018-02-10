volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_gh));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
);

MRF.makeRelative(phiHbyA);
adjustPhi(phiHbyA, U, p_gh);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_gh, U, phiHbyA, rAUf);

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
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_gh.relax();

p = p_gh + (g & mesh.C());

U = HbyA - rAU*fvc::grad(p_gh);
U.correctBoundaryConditions();
fvOptions.correct(U);
