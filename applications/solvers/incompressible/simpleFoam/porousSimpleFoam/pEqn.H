volVectorField HbyA("HbyA", U);
if (pressureImplicitPorosity)
{
    HbyA = trTU()&UEqn().H();
}
else
{
    HbyA = trAU()*UEqn().H();
}

UEqn.clear();
surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf());

mrfZones.makeRelative(phiHbyA);

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
