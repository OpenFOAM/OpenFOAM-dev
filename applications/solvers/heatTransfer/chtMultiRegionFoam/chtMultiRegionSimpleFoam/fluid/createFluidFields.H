    // Initialise fluid field pointer lists
    PtrList<rhoThermo> thermoFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<radiation::radiationModel> radiation(fluidRegions.size());

    List<scalar> initialMassFluid(fluidRegions.size());
    List<label> pRefCellFluid(fluidRegions.size(), 0);
    List<scalar> pRefValueFluid(fluidRegions.size(), 0.0);
    List<bool> frozenFlowFluid(fluidRegions.size(), false);

    PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
    PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

    PtrList<fv::IOoptionList> fluidFvOptions(fluidRegions.size());

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid\n" << endl;

        thermoFluid.set
        (
            i,
            rhoThermo::New(fluidRegions[i]).ptr()
        );

        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermoFluid[i].rho()
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*UFluid[i])
                    & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            compressible::turbulenceModel::New
            (
                rhoFluid[i],
                UFluid[i],
                phiFluid[i],
                thermoFluid[i]
            ).ptr()
        );

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
        );

        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        // Force p_rgh to be consistent with p
        p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

        radiation.set
        (
            i,
            radiation::radiationModel::New(thermoFluid[i].T())
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        const dictionary& simpleDict =
            fluidRegions[i].solutionDict().subDict("SIMPLE");

        setRefCell
        (
            thermoFluid[i].p(),
            p_rghFluid[i],
            simpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );

        simpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

        rhoMax.set(i, new dimensionedScalar(simpleDict.lookup("rhoMax")));
        rhoMin.set(i, new dimensionedScalar(simpleDict.lookup("rhoMin")));

        Info<< "    Adding fvOptions\n" << endl;
        fluidFvOptions.set
        (
            i,
            new fv::IOoptionList(fluidRegions[i])
        );
    }


