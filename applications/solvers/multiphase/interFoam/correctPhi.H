CorrectPhi
(
    U,
    phi,
    p_rgh,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
