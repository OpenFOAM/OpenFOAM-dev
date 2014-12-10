reducedUnits refUnits;

IOobject reducedUnitsDictIOobject
(
    "reducedUnitsDict",
    runTime.system(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
);

if (reducedUnitsDictIOobject.headerOk())
{
    Info<< nl
        << "Reading reference quantities from reducedUnitsDict file." << endl;

    IOdictionary reducedUnitsDict(reducedUnitsDictIOobject);

    refUnits.setRefValues(reducedUnitsDict);
}

Info << refUnits << endl;
