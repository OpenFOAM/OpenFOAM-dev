// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & Uf();

CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU(),
    pimple,
    true
);

#include "continuityErrs.H"

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
