Info<< "Reading CTT Properties\n" << endl;

// Your Alpha Constant
dimensionedScalar alphaAlpha
(
    "alpha",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), // Dimensionless
    0.0302011
);

// Velocity Field
volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// Pressure Field
volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

#include "createPhi.H"
