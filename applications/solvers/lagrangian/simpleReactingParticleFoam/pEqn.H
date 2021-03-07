// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
tUEqn.clear();
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p)
     ==
        clouds.Srho()
      + fvModels.source(psi, p, rho.name())
    );

    pEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pEqn.flux();
    }
}

p.relax();

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvConstraints.constrain(U);

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

Info<< "p min/max = " << min(p).value() << ", " << max(p).value() << endl;
