    const fvMesh& mesh = fluidRegions[i];

    CombustionModel<rhoReactionThermo>& reaction = reactionFluid[i];

    rhoReactionThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    label inertIndex = -1;
    if (Y.size())
    {
        const word inertSpecie(thermo.lookup("inertSpecie"));
        if (!composition.species().found(inertSpecie))
        {
            FatalIOErrorIn(args.executable().c_str(), thermo)
                << "Inert specie " << inertSpecie
                << " not found in available species "
                << composition.species()
                << exit(FatalIOError);
        }
        inertIndex = composition.species()[inertSpecie];
    }

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turbulence = turbulenceFluid[i];
    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsFluid[i];

    radiationModel& rad = radiation[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pimpleNoLoopControl& pimple = pimples.pimple(i);

    pressureControl& pressureControl = pressureControlFluid[i];

    scalar cumulativeContErr = cumulativeContErrs[i];
