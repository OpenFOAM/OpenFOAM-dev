CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
