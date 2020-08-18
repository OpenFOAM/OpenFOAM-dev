rho = thermo.rho();

const volScalarField rAU("rAU", 1.0/UEqn.A());
const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

tUEqn.clear();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvScalarMatrix p_rghEqn(p_rgh, dimMass/dimTime);

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p_rgh)*phiHbyA/fvc::interpolate(rho);

    while (simple.correctNonOrthogonal())
    {
        p_rghEqn =
            fvc::div(phiHbyA) + fvm::div(phid, p_rgh)
          - fvm::laplacian(rhorAUf, p_rgh)
         ==
            fvOptions(psi, p_rgh, rho.name());

        // Relax the pressure equation to ensure diagonal-dominance
        p_rghEqn.relax();

        p_rghEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        p_rghEqn.solve();
    }
}
else
{
    while (simple.correctNonOrthogonal())
    {
        p_rghEqn =
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p_rgh)
          ==
            fvOptions(psi, p_rgh, rho.name());

        p_rghEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        p_rghEqn.solve();
    }
}

phi = phiHbyA + p_rghEqn.flux();

p = p_rgh + rho*gh + pRef;

#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_rgh.relax();

// Correct the momentum source with the pressure gradient flux
// calculated from the relaxed pressure
U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
U.correctBoundaryConditions();
fvOptions.correct(U);

pressureControl.limit(p);

// For closed-volume compressible cases adjust the pressure level
// to obey overall mass continuity
if (closedVolume && !thermo.incompressible())
{
    p += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(psi);
    p_rgh = p - rho*gh - pRef;
    p.correctBoundaryConditions();
}

rho = thermo.rho();

if (!simple.transonic())
{
    rho.relax();
}
