CorrectPhi
(
    U,
    phi,
    p_gh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
