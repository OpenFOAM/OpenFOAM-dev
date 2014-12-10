Info<< "Reading thermal properties\n" << endl;

IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Switch thermalStress(thermalProperties.lookup("thermalStress"));

volScalarField threeKalpha
(
    IOobject
    (
        "threeKalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 0.0)
);


volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -1 , 0, 0), 0.0)
);


if (thermalStress)
{

    autoPtr<volScalarField> CPtr;

    IOobject CIO
    (
        "C",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& CDict(thermalProperties.subDict("C"));
    word CType(CDict.lookup("type"));
    if (CType == "uniform")
    {
        scalar CValue(readScalar(CDict.lookup("value")));

        CPtr.reset
        (
            new volScalarField
            (
                CIO,
                mesh,
                dimensionedScalar
                (
                    "C",
                    dimensionSet(0, 2, -2 , -1, 0),
                    CValue
                ),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );

    }
    else if (CType == "field")
    {
        CIO.readOpt() = IOobject::MUST_READ;

        CPtr.reset
        (
            new volScalarField
            (
                CIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorIn
        (
             "readThermalProperties.H"
        )   << "Valid type entries are uniform or field for C"
            << abort(FatalError);
    }

    volScalarField& C = CPtr();

    autoPtr<volScalarField> rhoKPtr;

    IOobject rhoKIO
    (
        "k",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& kDict(thermalProperties.subDict("k"));
    word kType(kDict.lookup("type"));
    if (kType == "uniform")
    {
        scalar rhoKValue(readScalar(kDict.lookup("value")));

        rhoKPtr.reset
        (
            new volScalarField
            (
                rhoKIO,
                mesh,
                dimensionedScalar
                (
                    "rhoK",
                    dimensionSet(1, 1, -3 , -1, 0),
                    rhoKValue
                ),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );

    }
    else if (kType == "field")
    {
        rhoKIO.readOpt() = IOobject::MUST_READ;

        rhoKPtr.reset
        (
            new volScalarField
            (
                rhoKIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorIn
        (
             "readThermalProperties.H"
        )   << "Valid type entries are uniform or field for K"
            << abort(FatalError);
    }

    volScalarField& rhoK = rhoKPtr();

    autoPtr<volScalarField> alphaPtr;

    IOobject alphaIO
    (
        "alpha",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );


    const dictionary& alphaDict(thermalProperties.subDict("alpha"));
    word alphaType(alphaDict.lookup("type"));

    if (alphaType == "uniform")
    {
        scalar alphaValue(readScalar(alphaDict.lookup("value")));
        alphaPtr.reset
        (
            new volScalarField
            (
                alphaIO,
                mesh,
                dimensionedScalar
                (
                    "alpha",
                    inv(dimTemperature),
                    alphaValue
                ),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
    }
    else if (alphaType == "field")
    {
        alphaIO.readOpt() = IOobject::MUST_READ;

        alphaPtr.reset
        (
            new volScalarField
            (
                alphaIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "readThermalProperties.H"
        )   << "Valid type entries are uniform or field for alpha"
            << abort(FatalError);
    }

    volScalarField& alpha = alphaPtr();

    Info<< "Normalising k : k/rho\n" << endl;
    volScalarField k(rhoK/rho);

    Info<< "Calculating thermal coefficients\n" << endl;

    threeKalpha = threeK*alpha;
    DT = k/C;
}
