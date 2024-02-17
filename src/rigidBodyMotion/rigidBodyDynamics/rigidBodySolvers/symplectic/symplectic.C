/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "symplectic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace rigidBodySolvers
{
    defineTypeNameAndDebug(symplectic, 0);
    addToRunTimeSelectionTable(rigidBodySolver, symplectic, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolvers::symplectic::symplectic
(
    rigidBodyMotion& body,
    const dictionary& dict
)
:
    rigidBodySolver(body)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolvers::symplectic::~symplectic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBD::rigidBodySolvers::symplectic::solve
(
    const scalarField& tau,
    const Field<spatialVector>& fx
)
{
    // First symplectic step:
    //     Half-step for linear and angular velocities
    //     Update position and orientation
    qDot() = qDot0() + 0.5*deltaT0()*qDdot();
    q() = q0() + deltaT()*qDot();

    correctQuaternionJoints();

    // Update the body-state prior to the evaluation of the restraints
    model_.forwardDynamicsCorrection(state());

    // Accumulate the restraint forces
    scalarField rtau(tau);
    Field<spatialVector> rfx(fx);
    model_.applyRestraints(rtau, rfx, state());

    // Calculate the body acceleration for the given state
    // and restraint forces
    model_.forwardDynamics(state(), rtau, rfx);

    // Second symplectic step:
    //     Complete update of linear and angular velocities
    qDot() += 0.5*deltaT()*qDdot();
}


// ************************************************************************* //
