    volScalarField rDeltaT
    (
        IOobject
        (
            "rDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless/dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    );
