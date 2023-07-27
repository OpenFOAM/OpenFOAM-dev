Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiuMulticomponentThermo> pThermo
(
    psiuMulticomponentThermo::New(mesh)
);
psiuMulticomponentThermo& thermo = pThermo();
thermo.validate(args.executable(), "ha", "ea");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.renameRho()
);

volScalarField& p = thermo.p();

volScalarField& b = thermo.Y("b");
Info<< "min(b) = " << min(b).value() << endl;

Info<< "\nReading field U\n" << endl;
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

#include "compressibleCreatePhi.H"

mesh.schemes().setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::RASModel> turbulence
(
    compressible::New<compressible::RASModel>
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
turbulenceThermophysicalTransportModels::unityLewisEddyDiffusivity
<
    RASThermophysicalTransportModel
    <
        ThermophysicalTransportModel
        <
            compressibleMomentumTransportModel,
            fluidThermo
        >
    >
> thermophysicalTransport(turbulence(), thermo, true);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.name(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


Info<< "Creating the unstrained laminar flame speed\n" << endl;
autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed
(
    laminarFlameSpeed::New(thermo)
);


Info<< "Reading strained laminar flame speed field Su\n" << endl;
volScalarField Su
(
    IOobject
    (
        "Su",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field betav\n" << endl;
volScalarField betav
(
    IOobject
    (
        "betav",
        mesh.facesInstance(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field Lobs\n" << endl;
volScalarField Lobs
(
    IOobject
    (
        "Lobs",
        mesh.facesInstance(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field CT\n" << endl;
volSymmTensorField CT
(
    IOobject
    (
        "CT",
        mesh.facesInstance(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field Nv\n" << endl;
volScalarField Nv
(
    IOobject
    (
        "Nv",
        mesh.facesInstance(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field nsv\n" << endl;
volSymmTensorField nsv
(
    IOobject
    (
        "nsv",
        mesh.facesInstance(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

IOdictionary PDRProperties
(
    IOobject
    (
        "PDRProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//- Create the drag model
autoPtr<PDRDragModel> drag = PDRDragModel::New
(
    PDRProperties,
    turbulence,
    rho,
    U,
    phi
);

//- Create the flame-wrinkling model
autoPtr<XiModel> flameWrinkling = XiModel::New
(
    PDRProperties,
    thermo,
    turbulence,
    Su,
    rho,
    b,
    phi
);

Info<< "Calculating turbulent flame speed field St\n" << endl;
volScalarField St
(
    IOobject
    (
        "St",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    flameWrinkling->Xi()*Su
);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

if (thermo.containsSpecie("ft"))
{
    fields.add(thermo.Y("ft"));
}

fields.add(b);
fields.add(thermo.he());
fields.add(thermo.heu());
flameWrinkling->addXi(fields);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
