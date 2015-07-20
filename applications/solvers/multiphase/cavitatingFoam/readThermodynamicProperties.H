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
    thermodynamicProperties.lookup("psil")
);

dimensionedScalar rholSat
(
    "rholSat",
    dimDensity,
    thermodynamicProperties.lookup("rholSat")
);

dimensionedScalar psiv
(
    "psiv",
    dimCompressibility,
    thermodynamicProperties.lookup("psiv")
);

dimensionedScalar pSat
(
    "pSat",
    dimPressure,
    thermodynamicProperties.lookup("pSat")
);

dimensionedScalar rhovSat("rhovSat", psiv*pSat);

dimensionedScalar rhol0("rhol0", rholSat - pSat*psil);

dimensionedScalar rhoMin
(
    "rhoMin",
    dimDensity,
    thermodynamicProperties.lookup("rhoMin")
);
