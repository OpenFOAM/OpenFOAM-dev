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

if (fluid.found("pMin"))
{
    IOWarningInFunction(fluid)
        << "Pressure limits, pMin and pMax, are now read from "
        << pimple.dict().name() << endl;
}

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    fluid.incompressible()
);

if (fluid.incompressible())
{
    p = p_rgh + fluid.rho()*gh;
}

if (p_rgh.needReference() && fluid.incompressible())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
}

p_rgh = p - fluid.rho()*gh;

mesh.schemes().setFluxRequired(p_rgh.name());

PtrList<volScalarField> rAUs;
PtrList<surfaceScalarField> rAUfs;
