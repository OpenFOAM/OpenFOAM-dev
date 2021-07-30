Info<< "Reading field p (kinematic)\n" << endl;
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


autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

#include "readGravitationalAcceleration.H"

Info<< "Creating field zeta\n" << endl;
volVectorField zeta
(
    IOobject
    (
        "zeta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength, Zero)
);

Info<< "Creating field p_gh\n" << endl;
volScalarField p_gh
(
    IOobject
    (
        "p_gh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_gh to be consistent with p
// Height is made relative to field 'refLevel'
p_gh =  p - (g & mesh.C());

pressureReference pressureReference(p_gh, pimple.dict());

mesh.setFluxRequired(p_gh.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
