tmp<volScalarField> rAU;

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime/dimDensity, 1)
    );

    correctUphiBCs(U, phi, true);

    CorrectPhi
    (
        phi,
        U,
        p_rgh,
        surfaceScalarField("rAUf", fvc::interpolate(rAU())),
        geometricZeroField(),
        pressureReference,
        pimple
    );

    #include "continuityErrs.H"
}
else
{
    correctUphiBCs(U, phi, true);

    CorrectPhi
    (
        phi,
        U,
        p_rgh,
        dimensionedScalar(dimTime/rho.dimensions(), 1),
        geometricZeroField(),
        pressureReference,
        pimple
    );

    #include "continuityErrs.H"
}
