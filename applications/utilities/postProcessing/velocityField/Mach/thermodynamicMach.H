#   include "readThermodynamicProperties.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Theader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );


        // Check U exists
        if (Uheader.headerOk() && Theader.headerOk())
        {
            mesh.readUpdate();

            volVectorField U(Uheader, mesh);
            volScalarField T(Theader, mesh);

            volScalarField Ma
            (
                IOobject
                (
                    "Ma",
                    runTime.timeName(),
                    mesh
                ),
                mag(U)/(sqrt(((Cv + R)/Cv)*R*T))
            );
            Ma.write();
        }
        else
        {
            Info<< "    No U or T" << endl;
        }
    }
