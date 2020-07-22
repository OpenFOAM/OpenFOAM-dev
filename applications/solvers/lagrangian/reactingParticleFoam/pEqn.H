if (!pimple.simpleRho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::interpolate(rho)*fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
    )
);

fvc::makeRelative(phiHbyA, rho, U);
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

fvScalarMatrix pDDtEqn
(
    fvc::ddt(rho) + psi*correction(fvm::ddt(p))
  + fvc::div(phiHbyA)
 ==
    clouds.Srho()
  + surfaceFilm.Srho()
  + fvOptions(psi, p, rho.name())
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        pDDtEqn
      - fvm::laplacian(rhorAUf, p)
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pEqn.flux();
    }
}

p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

bool limitedp = pressureControl.limit(p);

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (limitedp)
{
    rho = thermo.rho();
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

if (pimple.simpleRho())
{
    rho = thermo.rho();
}

// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, rho, U, phi);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
