/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "singleRigidBodyMeshMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleRigidBodyMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        singleRigidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::List<Foam::septernion>
Foam::singleRigidBodyMeshMotion::transforms0() const
{
    return List<septernion>(1, SBMFPtr_().transformation());
}


void Foam::singleRigidBodyMeshMotion::moveBodies()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleRigidBodyMeshMotion::singleRigidBodyMeshMotion
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    multiRigidBodyMeshMotion(name, mesh, dict),
    SBMFPtr_(solidBodyMotionFunction::New(dict, mesh.time(), "function"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleRigidBodyMeshMotion::~singleRigidBodyMeshMotion()
{}


// ************************************************************************* //
