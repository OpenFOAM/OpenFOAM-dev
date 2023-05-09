/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solid.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMeshMover.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(solid, 0);
    addToRunTimeSelectionTable(solver, solid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::solid::readControls()
{
    maxDi =
        runTime.controlDict().lookupOrDefault<scalar>("maxDi", 1.0);

    maxDeltaT_ =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", vGreat);
}


void Foam::solvers::solid::correctDiNum()
{
    const volScalarField kappa
    (
        thermo.isotropic()
      ? thermo.kappa()
      : mag(thermo.Kappa())()
    );

    const volScalarField::Internal DiNumvf
    (
        fvc::surfaceSum
        (
            mesh.magSf()
           *fvc::interpolate(kappa)
           *mesh.surfaceInterpolation::deltaCoeffs()
        )()()
       /(mesh.V()*thermo.rho()()*thermo.Cp()())
       *runTime.deltaT()
    );

    const scalar meanDiNum = gAverage(DiNumvf);
    const scalar maxDiNum = gMax(DiNumvf);

    Info<< "Diffusion Number mean: " << meanDiNum
        << " max: " << maxDiNum << endl;

    DiNum = maxDiNum;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::solid::solid
(
    fvMesh& mesh,
    autoPtr<solidThermo> thermoPtr
)
:
    solver(mesh),

    thermoPtr_(thermoPtr),
    thermo_(thermoPtr_()),

    T_(thermo_.T()),

    thermophysicalTransport(solidThermophysicalTransportModel::New(thermo_)),

    DiNum(0),

    thermo(thermo_),
    T(T_)
{
    thermo.validate("solid", "h", "e");

    if (transient())
    {
        correctDiNum();
    }
    else if (LTS)
    {
        FatalError
            << type()
            << " solver does not support LTS, use 'steadyState' ddtScheme"
            << exit(FatalError);
    }
}



Foam::solvers::solid::solid(fvMesh& mesh)
:
    solid(mesh, solidThermo::New(mesh))
{
    // Read the controls
    readControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::solid::~solid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::solid::maxDeltaT() const
{
    scalar deltaT = min(fvModels().maxDeltaT(), maxDeltaT_);

    if (DiNum > small)
    {
        deltaT = min(deltaT, maxDi/DiNum*runTime.deltaTValue());
    }

    return deltaT;
}


void Foam::solvers::solid::preSolve()
{
    // Read the controls
    readControls();

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();

    if (transient())
    {
        correctDiNum();
    }
}


void Foam::solvers::solid::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        if (!mesh_.mover().solidBody())
        {
            FatalErrorInFunction
                << "Region " << name() << " of type " << type()
                << " does not support non-solid body mesh motion"
                << exit(FatalError);
        }

        mesh_.move();
    }
}


void Foam::solvers::solid::prePredictor()
{
    if (pimple.predictTransport())
    {
        thermophysicalTransport->predict();
    }
}


void Foam::solvers::solid::momentumPredictor()
{}


void Foam::solvers::solid::thermophysicalPredictor()
{
    volScalarField& e = thermo_.he();
    const volScalarField& rho = thermo.rho();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho, e)
          + thermophysicalTransport->divq(e)
          ==
            fvModels().source(rho, e)
        );

        eEqn.relax();

        fvConstraints().constrain(eEqn);

        eEqn.solve();

        fvConstraints().constrain(e);
    }

    thermo_.correct();
}


void Foam::solvers::solid::pressureCorrector()
{}


void Foam::solvers::solid::postCorrector()
{
    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
}


void Foam::solvers::solid::postSolve()
{}


// ************************************************************************* //
