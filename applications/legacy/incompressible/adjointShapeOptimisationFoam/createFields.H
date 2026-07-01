Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
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
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


pressureReference pressureReference(p, simple.dict());

mesh.schemes().setFluxRequired(p.name());


Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhia.H"


label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.schemes().setFluxRequired(pa.name());


autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);


const dimensionedScalar zeroSensitivity("0", sqr(dimensions::velocity), 0.0);
const dimensionedScalar zeroAlpha("0", dimensions::rate, 0.0);

dimensionedScalar lambda
(
    "lambda",
    inv(dimensions::kinematicViscosity),
    viscosity
);
const dimensionedScalar alphaMax("alphaMax", dimensions::rate, viscosity);

const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
//const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    lambda*max(Ua & U, zeroSensitivity)
);
zeroCells(alpha, inletCells);
//zeroCells(alpha, outletCells);

#include "createFvModels.H"
#include "createFvConstraints.H"
