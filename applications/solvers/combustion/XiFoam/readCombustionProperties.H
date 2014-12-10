    Info<< "Reading combustion properties\n" << endl;

    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    word SuModel
    (
        combustionProperties.lookup("SuModel")
    );

    dimensionedScalar sigmaExt
    (
        combustionProperties.lookup("sigmaExt")
    );

    word XiModel
    (
        combustionProperties.lookup("XiModel")
    );

    dimensionedScalar XiCoef
    (
        combustionProperties.lookup("XiCoef")
    );

    dimensionedScalar XiShapeCoef
    (
        combustionProperties.lookup("XiShapeCoef")
    );

    dimensionedScalar uPrimeCoef
    (
        combustionProperties.lookup("uPrimeCoef")
    );

    ignition ign(combustionProperties, runTime, mesh);
