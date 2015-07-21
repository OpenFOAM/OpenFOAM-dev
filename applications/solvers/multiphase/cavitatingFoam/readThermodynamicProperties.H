Info<< "Reading thermodynamicProperties\n" << endl;

IOdictionary thermodynamicProperties
(
    IOobject
    (
        "thermodynamicProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar psil
(
    "psil",
    dimCompressibility,
    thermodynamicProperties
);

dimensionedScalar rholSat
(
    "rholSat",
    dimDensity,
    thermodynamicProperties
);

dimensionedScalar psiv
(
    "psiv",
    dimCompressibility,
    thermodynamicProperties
);

dimensionedScalar pSat
(
    "pSat",
    dimPressure,
    thermodynamicProperties
);

dimensionedScalar rhovSat("rhovSat", psiv*pSat);

dimensionedScalar rhol0("rhol0", rholSat - pSat*psil);

dimensionedScalar rhoMin
(
    "rhoMin",
    dimDensity,
    thermodynamicProperties
);
