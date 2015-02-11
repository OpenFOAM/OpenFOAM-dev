CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    divU,
    pimple
);

//***HGW phi.oldTime() = phi;

#include "continuityErrs.H"
