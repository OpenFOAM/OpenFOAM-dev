Info<< "Reading financial properties\n" << endl;

IOdictionary financialProperties
(
    IOobject
    (
        "financialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar strike
(
    "strike",
    dimensions::length,
    financialProperties
);
const dimensionedScalar r("r", dimless/dimensions::time, financialProperties);

const dimensionedScalar sigma
(
    "sigma",
    dimensionSet(0, 0, -0.5, 0, 0),
    financialProperties
);

const dimensionedScalar sigmaSqr = sqr(sigma);


Info<< nl << "Reading field V" << endl;

volScalarField V
(
    IOobject
    (
        "V",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


surfaceVectorField Pf
(
    IOobject
    (
        "Pf",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.Cf()
);


volVectorField P
(
    IOobject
    (
        "P",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.C()
);

V == max
(
    P.component(Foam::vector::X) - strike,
    dimensionedScalar(V.dimensions(), 0)
);

volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(V)().component(Foam::vector::X)
);
