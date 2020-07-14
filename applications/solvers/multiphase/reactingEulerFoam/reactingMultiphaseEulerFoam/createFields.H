#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<phaseSystem> fluidPtr
(
    phaseSystem::New(mesh)
);
phaseSystem& fluid = fluidPtr();
phaseSystem::phaseModelList& phases = fluid.phases();

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
if (fluid.incompressible())
{
    p = max(p_rgh + fluid.rho()*gh, pMin);

    if (p_rgh.needReference())
    {
        setRefCell
        (
            p,
            p_rgh,
            pimple.dict(),
            pRefCell,
            pRefValue
        );

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - fluid.rho()*gh;
    }
}
mesh.setFluxRequired(p_rgh.name());

PtrList<volScalarField> rAUs;
PtrList<surfaceScalarField> rAUfs;
