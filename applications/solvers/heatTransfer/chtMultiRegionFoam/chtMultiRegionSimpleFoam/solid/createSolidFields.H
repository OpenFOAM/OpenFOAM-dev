// Initialise solid field pointer lists
PtrList<coordinateSystem> coordinates(solidRegions.size());
PtrList<solidThermo> thermos(solidRegions.size());
PtrList<radiation::radiationModel> radiations(solidRegions.size());
PtrList<fv::options> solidHeatSources(solidRegions.size());
PtrList<volScalarField> betavSolid(solidRegions.size());
PtrList<volSymmTensorField> aniAlphas(solidRegions.size());

List<bool> residualReachedSolid(solidRegions.size(), true);
List<bool> residualControlUsedSolid(solidRegions.size(), false);
List<bool> firstIterationSolid(solidRegions.size(), true);

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
        new fv::options(solidRegions[i])
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
                    Zero
                ),
                zeroGradientFvPatchSymmTensorField::typeName
            )
        );

        aniAlphas[i].primitiveFieldRef() =
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

    if (betavSolidIO.typeHeaderOk<volScalarField>(true))
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
                dimensionedScalar("1", dimless, scalar(1))
            )
        );
    }

    const dictionary& simpleDict =
        solidRegions[i].solutionDict().subDict("SIMPLE");

    if (simpleDict.isDict("residualControl"))
    {
        Info<< "    Reading residual controls\n" << endl;

        bool& residualControlUsed = residualControlUsedSolid[i];
        bool& resReachedSolid = residualReachedSolid[i];
        bool& firstIteration = firstIterationSolid[i];
        const fvMesh& mesh = solidRegions[i];

        #include "readSolidMultiRegionResidualControls.H"
    }
}
