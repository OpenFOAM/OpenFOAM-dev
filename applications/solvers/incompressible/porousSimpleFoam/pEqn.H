tmp<volVectorField> tHbyA;
if (pressureImplicitPorosity)
{
    tHbyA = constrainHbyA(trTU()&UEqn.H(), U, p);
}
else
{
    tHbyA = constrainHbyA(trAU()*UEqn.H(), U, p);
}
volVectorField& HbyA = tHbyA.ref();

tUEqn.clear();
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

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

    fvScalarMatrix pEqn = tpEqn.ref();

    pEqn.setReference
    (
        pressureReference.refCell(),
        pressureReference.refValue()
    );

    pEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
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
fvConstraints.constrain(U);
