#include "readThermophysicalProperties.H"

Info<< "Reading displacement field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/setting displacement correction field Dcorr\n" << endl;
volVectorField Dcorr
(
    IOobject
    (
        "Dcorr",
        runTime.name(),
        mesh
    ),
    D
);
Dcorr = Zero;

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(fvc::grad(D)) + (lambda*I)*tr(fvc::grad(D))
);

volTensorField sigmaExp
(
    IOobject
    (
        "sigmaExp",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (lambda - mu)*fvc::grad(Dcorr) + mu*fvc::grad(Dcorr)().T()
  + (lambda*I)*tr(fvc::grad(Dcorr))
);
