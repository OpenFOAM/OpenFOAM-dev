    Info<< "\nReading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    word polyatomicCloudName("polyatomicCloud");

    potential polyPot
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                polyatomicCloudName +  "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    polyatomicCloud polyatomics(polyatomicCloudName, mesh, polyPot);

    word monoatomicCloudName("monoatomicCloud");

    potential monoPot
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                monoatomicCloudName +  "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    monoatomicCloud monoatomics(monoatomicCloudName, mesh, monoPot);
