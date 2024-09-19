/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::solid::dependenciesModified() const
{
    return runTime.controlDict().modified();
}


bool Foam::solvers::solid::read()
{
    solver::read();

    maxDeltaT_ =
        runTime.controlDict().found("maxDeltaT")
      ? runTime.controlDict().lookup<scalar>("maxDeltaT", runTime.userUnits())
      : vGreat;

    return true;
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

    thermo(thermo_),
    T(T_)
{
    thermo.validate("solid", "h", "e");

    if (LTS)
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
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::solid::~solid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::solid::maxDeltaT() const
{
    return min(fvModels().maxDeltaT(), maxDeltaT_);
}


void Foam::solvers::solid::preSolve()
{
    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
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


void Foam::solvers::solid::motionCorrector()
{}


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
    const volScalarField& rho = thermo_.rho();

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

        thermo_.correct();
    }
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
