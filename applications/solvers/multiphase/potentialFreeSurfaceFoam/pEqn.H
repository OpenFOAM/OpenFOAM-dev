rAU = 1.0/UEqn.A();
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_gh));

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi, Uf))
);

MRF.makeRelative(phiHbyA);

if (p_gh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p_gh);
    fvc::makeAbsolute(phiHbyA, U);
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_gh, U, phiHbyA, rAUf);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_ghEqn
    (
        fvm::laplacian(rAUf, p_gh) == fvc::div(phiHbyA)
    );

    p_ghEqn.setReference
    (
        pressureReference.refCell(),
        pressureReference.refValue()
    );

    p_ghEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_ghEqn.flux();

        // Explicitly relax pressure for momentum corrector
        p_gh.relax();

        U = HbyA - rAU*fvc::grad(p_gh);
        U.correctBoundaryConditions();
        fvConstraints.constrain(U);
    }
}

#include "continuityErrs.H"

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

p = p_gh + (g & mesh.C());
