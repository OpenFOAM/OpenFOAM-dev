#include "readThermodynamicProperties.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
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
);

#include "createPhi.H"

mesh.schemes().setFluxRequired(p.name());

// Mass flux (corrected by rhoEqn.H)
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

Info<< "Reading phaseProperties\n" << endl;

incompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alphav(mixture.alpha1());
alphav.oldTime();

volScalarField& alphal(mixture.alpha2());

Info<< "Creating compressibilityModel\n" << endl;
autoPtr<barotropicCompressibilityModel> psiModel =
    barotropicCompressibilityModel::New
    (
        thermodynamicProperties,
        alphav
    );

const volScalarField& psi = psiModel->psi();

rho == max
(
    psi*p
  + alphal*rhol0
  + ((alphav*psiv + alphal*psil) - psi)*pSat,
    rhoMin
);

mesh.schemes().setFluxRequired(p.name());
mesh.schemes().setFluxRequired(rho.name());

// Create incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);
