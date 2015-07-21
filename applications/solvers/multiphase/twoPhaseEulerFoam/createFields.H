#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating twoPhaseSystem\n" << endl;

twoPhaseSystem fluid(mesh, g);

phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();

volScalarField& alpha1 = phase1;
volScalarField& alpha2 = phase2;

volVectorField& U1 = phase1.U();
surfaceScalarField& phi1 = phase1.phi();
surfaceScalarField& alphaPhi1 = phase1.alphaPhi();
surfaceScalarField& alphaRhoPhi1 = phase1.alphaRhoPhi();

volVectorField& U2 = phase2.U();
surfaceScalarField& phi2 = phase2.phi();
surfaceScalarField& alphaPhi2 = phase2.alphaPhi();
surfaceScalarField& alphaRhoPhi2 = phase2.alphaRhoPhi();

surfaceScalarField& phi = fluid.phi();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

rhoThermo& thermo1 = phase1.thermo();
rhoThermo& thermo2 = phase2.thermo();

volScalarField& p = thermo1.p();

volScalarField& rho1 = thermo1.rho();
const volScalarField& psi1 = thermo1.psi();

volScalarField& rho2 = thermo2.rho();
const volScalarField& psi2 = thermo2.psi();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);


Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));
