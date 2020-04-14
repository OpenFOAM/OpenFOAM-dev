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


Info<< "Creating face flux\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(mesh.Sf().dimensions()*U.dimensions(), 0)
);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

dimensionedVector Ubar("Ubar", dimVelocity, laminarTransport);

vector flowDirection = (Ubar/mag(Ubar)).value();
tensor flowMask = sqr(flowDirection);

dimensionedVector gradP
(
    "gradP",
    dimensionSet(0, 1, -2, 0, 0),
    Zero
);

#include "createFvOptions.H"
