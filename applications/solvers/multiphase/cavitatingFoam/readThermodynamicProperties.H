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

    dimensionedScalar psil(thermodynamicProperties.lookup("psil"));

    dimensionedScalar rholSat(thermodynamicProperties.lookup("rholSat"));

    dimensionedScalar psiv(thermodynamicProperties.lookup("psiv"));

    dimensionedScalar pSat(thermodynamicProperties.lookup("pSat"));

    dimensionedScalar rhovSat("rhovSat", psiv*pSat);

    dimensionedScalar rhol0("rhol0", rholSat - pSat*psil);

    dimensionedScalar rhoMin(thermodynamicProperties.lookup("rhoMin"));
