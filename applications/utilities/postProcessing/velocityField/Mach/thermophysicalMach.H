    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (Uheader.headerOk())
        {
            volVectorField U(Uheader, mesh);

            autoPtr<fluidThermo> thermo
            (
                fluidThermo::New(mesh)
            );

            volScalarField Cp = thermo->Cp();
            volScalarField Cv = thermo->Cv();

            volScalarField Ma
            (
                IOobject
                (
                    "Ma",
                    runTime.timeName(),
                    mesh
                ),
                mag(U)/(sqrt((Cp/Cv)*(Cp - Cv)*thermo->T()))
            );
            Ma.write();
        }
        else
        {
            Info<< "    No U" << endl;
        }
    }
