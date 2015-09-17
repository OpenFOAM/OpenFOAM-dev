#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<multiphaseSystem> fluidPtr
(
    multiphaseSystem::New(mesh)
);
multiphaseSystem& fluid = fluidPtr();
multiphaseSystem::phaseModelList& phases = fluid.phases();

surfaceScalarField& phi = fluid.phi();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = phases[0].thermo().p();

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

const IOMRFZoneList& MRF = fluid.MRF();
fv::IOoptionList& fvOptions = fluid.fvOptions();
