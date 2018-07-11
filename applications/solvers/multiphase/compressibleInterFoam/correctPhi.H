CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    divU,
    pimple
);

#include "continuityErrs.H"
