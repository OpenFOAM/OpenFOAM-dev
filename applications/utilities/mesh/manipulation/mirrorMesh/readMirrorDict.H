    Info<< "Reading mirrorMeshDict\n" << endl;

    IOdictionary mirrorMeshDict
    (
        IOobject
        (
            "mirrorMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    plane mirrorPlane(mirrorMeshDict);

    scalar planeTolerance
    (
        readScalar(mirrorMeshDict.lookup("planeTolerance"))
    );
