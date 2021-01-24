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
