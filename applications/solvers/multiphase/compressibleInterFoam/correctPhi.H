// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & Uf();

correctUphiBCs(U, phi, true);

CorrectPhi
(
    phi,
    U,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU(),
    pimple
);

#include "continuityErrs.H"

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
