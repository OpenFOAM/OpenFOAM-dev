tmp<volScalarField> rAU;

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime/dimDensity, 1)
    );

    #include "correctPhi.H"
}
else
{
    CorrectPhi
    (
        U,
        phi,
        p_rgh,
        dimensionedScalar(dimTime/rho.dimensions(), 1),
        geometricZeroField(),
        pimple
    );

    #include "continuityErrs.H"
}
