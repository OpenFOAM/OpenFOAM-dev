#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<multiphaseSystem> fluidPtr
(
    multiphaseSystem::New(mesh)
);
multiphaseSystem& fluid = fluidPtr();
multiphaseSystem::phaseModelList& phases = fluid.phases();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = phases[0].thermoRef().p();

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
