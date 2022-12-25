Info<< "Reading phaseProperties\n" << endl;

cavitatingTwoPhaseMixture mixture(mesh);

volScalarField& p(mixture.p());
volScalarField& rho(mixture.rho());

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
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
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

mesh.schemes().setFluxRequired(p.name());
mesh.schemes().setFluxRequired(rho.name());

// Create incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);
