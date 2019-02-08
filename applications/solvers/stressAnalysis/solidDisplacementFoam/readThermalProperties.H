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
    dimensionedScalar(dimensionSet(0, 2, -2 , -1, 0), 0)
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
    dimensionedScalar(dimensionSet(0, 2, -1 , 0, 0), 0)
);


// Cache "k" for use in thermal BCs
autoPtr<volScalarField> rhoKPtr;

if (thermalStress)
{
    IOobject CIO
    (
        "C",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    autoPtr<volScalarField> CPtr;

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
                    dimensionSet(0, 2, -2 , -1, 0),
                    CValue
                )
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
        FatalErrorInFunction
            << "Valid type entries are uniform or field for C"
            << abort(FatalError);
    }

    volScalarField& C = CPtr();

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
                    dimensionSet(1, 1, -3 , -1, 0),
                    rhoKValue
                )
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
        FatalErrorInFunction
            << "Valid type entries are uniform or field for K"
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
                    inv(dimTemperature),
                    alphaValue
                )
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
        FatalErrorInFunction
            << "Valid type entries are uniform or field for alpha"
            << abort(FatalError);
    }

    volScalarField& alpha = alphaPtr();

    Info<< "Normalising k : k/rho\n" << endl;
    volScalarField k(rhoK/rho);

    Info<< "Calculating thermal coefficients\n" << endl;

    threeKalpha = threeK*alpha;
    DT = k/C;
}
