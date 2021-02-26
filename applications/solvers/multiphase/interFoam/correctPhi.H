// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & Uf();

if (divU.valid())
{
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
}
else
{
    CorrectPhi
    (
        U,
        phi,
        p_rgh,
        surfaceScalarField("rAUf", fvc::interpolate(rAU())),
        geometricZeroField(),
        pimple,
        true
    );
}

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phi, U);
