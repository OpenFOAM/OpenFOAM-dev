CorrectPhi
(
    U,
    phi,
    p,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU(),
    pimple,
    true
);

#include "continuityErrs.H"
