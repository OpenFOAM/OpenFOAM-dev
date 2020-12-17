    const fvMesh& mesh = fluidRegions[i];

    combustionModel& reaction = reactionFluid[i];

    fluidReactionThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::momentumTransportModel& turbulence = turbulenceFluid[i];
    fluidReactionThermophysicalTransportModel& thermophysicalTransport =
        thermophysicalTransportFluid[i];

    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const dimensionedScalar& pRef = pRefFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    const fv::options& fvOptions = fvOptionsFluid[i];

    #include "checkRadiationModel.H"

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pimpleNoLoopControl& pimple = pimples.pimple(i);

    pressureControl& pressureControl = pressureControlFluid[i];

    scalar cumulativeContErr = cumulativeContErrs[i];

    // This solver does not support moving mesh but it uses the pressure
    // equation of one which does, so we need a dummy face-momentum field
    autoPtr<surfaceVectorField> rhoUf(nullptr);
