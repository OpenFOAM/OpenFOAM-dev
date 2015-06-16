scalar maxDeltaT
(
    pimple.dict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
);

volScalarField rDeltaT
(
    IOobject
    (
        "rDeltaT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rSubDeltaT
(
    IOobject
    (
        "rSubDeltaT",
        runTime.timeName(),
        mesh
    ),
    mesh,
    1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT)
);
