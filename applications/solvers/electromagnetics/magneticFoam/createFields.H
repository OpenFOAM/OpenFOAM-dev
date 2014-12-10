    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    List<magnet> magnets(transportProperties.lookup("magnets"));

    surfaceScalarField murf
    (
        IOobject
        (
            "murf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        1
    );

    surfaceScalarField Mrf
    (
        IOobject
        (
            "Mrf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("Mr", dimensionSet(0, 1, 0, 0, 0, 1, 0), 0)
    );

    forAll(magnets, i)
    {
        label magnetZonei = mesh.faceZones().findZoneID(magnets[i].name());

        if (magnetZonei == -1)
        {
            FatalIOErrorIn(args.executable().c_str(), transportProperties)
                << "Cannot find faceZone for magnet " << magnets[i].name()
                << exit(FatalIOError);
        }

        const labelList& faces = mesh.faceZones()[magnetZonei];

        const scalar muri = magnets[i].mur();
        const scalar Mri = magnets[i].Mr().value();
        const vector& orientationi = magnets[i].orientation();

        const surfaceVectorField& Sf = mesh.Sf();

        forAll(faces, i)
        {
            label facei = faces[i];
            murf[facei] = muri;
            Mrf[facei] = Mri*(orientationi & Sf[facei]);
        }
    }
