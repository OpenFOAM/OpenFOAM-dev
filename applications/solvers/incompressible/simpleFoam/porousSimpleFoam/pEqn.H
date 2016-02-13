tmp<volVectorField> tHbyA;
if (pressureImplicitPorosity)
{
    tHbyA = constrainHbyA(trTU()&UEqn().H(), U, p);
}
else
{
    tHbyA = constrainHbyA(trAU()*UEqn().H(), U, p);
}
volVectorField& HbyA = tHbyA();

UEqn.clear();
surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf());

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn = (fvm::laplacian(trTU(), p) == fvc::div(phiHbyA));
    }
    else
    {
        tpEqn = (fvm::laplacian(trAU(), p) == fvc::div(phiHbyA));
    }

    tpEqn().setReference(pRefCell, pRefValue);

    tpEqn().solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - tpEqn().flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U = HbyA - (trTU()&fvc::grad(p));
}
else
{
    U = HbyA - (trAU()*fvc::grad(p));
}

U.correctBoundaryConditions();
fvOptions.correct(U);
