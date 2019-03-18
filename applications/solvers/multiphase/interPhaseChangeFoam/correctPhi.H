CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU,
    pimple,
    true
);

#include "continuityErrs.H"
