volScalarField rAUrel(1.0/UrelEqn.A());
volVectorField HbyA("HbyA", Urel);
HbyA = rAUrel*UrelEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + fvc::interpolate(rAUrel)*fvc::ddtCorr(Urel, phi)
);

adjustPhi(phiHbyA, Urel, p);

tmp<volScalarField> rAtUrel(rAUrel);

if (pimple.consistent())
{
    rAtUrel = 1.0/max(1.0/rAUrel - UrelEqn.H1(), 0.1/rAUrel);
    phiHbyA +=
        fvc::interpolate(rAtUrel() - rAUrel)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAUrel - rAtUrel())*fvc::grad(p);
}

if (pimple.nCorrPiso() <= 1)
{
    tUrelEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, Urel, phiHbyA, rAtUrel());

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtUrel(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

p.relax();

// Momentum corrector
Urel = HbyA - rAtUrel()*fvc::grad(p);
Urel.correctBoundaryConditions();
fvConstraints.constrain(Urel);
