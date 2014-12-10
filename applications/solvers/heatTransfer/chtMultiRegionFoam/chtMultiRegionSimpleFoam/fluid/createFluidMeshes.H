    const wordList fluidNames(rp["fluid"]);

    PtrList<fvMesh> fluidRegions(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info<< "Create fluid mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    fluidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
