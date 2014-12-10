Info<< nl << "Reading MD Equilibration Dictionary" << nl << endl;

IOdictionary mdEquilibrationDict
(
    IOobject
    (
        "mdEquilibrationDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar targetTemperature = readScalar
(
    mdEquilibrationDict.lookup("targetTemperature")
);
