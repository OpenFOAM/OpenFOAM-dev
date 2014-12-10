autoPtr<basicThermo> thermo
(
    basicThermo::New(mesh)
);

const volScalarField& h = thermo->he();

// Register copy of thermo density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo->rho()
);

// Construct turbulence model (if fluid)
autoPtr<volVectorField> UPtr;
autoPtr<surfaceScalarField> phiPtr;
autoPtr<compressible::turbulenceModel> turbulence;

if (isA<fluidThermo>(thermo()))
{
    UPtr.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    const volVectorField& U = UPtr();

    #include "compressibleCreatePhi.H"
    // Copy phi to autoPtr. Rename to make sure copy is now registered as 'phi'.
    phi.rename("phiFluid");
    phiPtr.reset(new surfaceScalarField("phi", phi));

    turbulence = compressible::turbulenceModel::New
    (
        rho,
        U,
        phiPtr(),
        refCast<const fluidThermo>(thermo())
    );
}
