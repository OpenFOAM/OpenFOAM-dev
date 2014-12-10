    // Initialise solid field pointer lists
    PtrList<coordinateSystem> coordinates(solidRegions.size());
    PtrList<solidThermo> thermos(solidRegions.size());
    PtrList<radiation::radiationModel> radiations(solidRegions.size());
    PtrList<fv::IOoptionList> solidHeatSources(solidRegions.size());
    PtrList<volScalarField> betavSolid(solidRegions.size());
    PtrList<volSymmTensorField> aniAlphas(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermos\n" << endl;
        thermos.set(i, solidThermo::New(solidRegions[i]));

        Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::IOoptionList(solidRegions[i])
        );

        if (!thermos[i].isotropic())
        {
            Info<< "    Adding coordinateSystems\n" << endl;
            coordinates.set
            (
                i,
                coordinateSystem::New(solidRegions[i], thermos[i])
            );

            tmp<volVectorField> tkappaByCp =
                thermos[i].Kappa()/thermos[i].Cp();

            aniAlphas.set
            (
                i,
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedSymmTensor
                    (
                        "zero",
                        tkappaByCp().dimensions(),
                        symmTensor::zero
                    ),
                    zeroGradientFvPatchSymmTensorField::typeName
                )
            );

            aniAlphas[i].internalField() =
                coordinates[i].R().transformVector(tkappaByCp());
            aniAlphas[i].correctBoundaryConditions();

        }

        IOobject betavSolidIO
        (
            "betavSolid",
            runTime.timeName(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (betavSolidIO.headerOk())
        {
            betavSolid.set
            (
                i,
                new volScalarField(betavSolidIO, solidRegions[i])
            );
        }
        else
        {
            betavSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("1", dimless, scalar(1.0))
                )
            );
        }
    }
