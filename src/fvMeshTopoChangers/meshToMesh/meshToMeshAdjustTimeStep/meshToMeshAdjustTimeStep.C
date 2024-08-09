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

#include "meshToMeshAdjustTimeStep.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(meshToMeshAdjustTimeStep, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        meshToMeshAdjustTimeStep,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::meshToMeshAdjustTimeStep::meshToMeshAdjustTimeStep
(
    const word& name,
    const objectRegistry& obr
)
:
    fvMeshFunctionObject(name, obr),
    meshToMesh_
    (
        refCast<const fvMeshTopoChangers::meshToMesh>(mesh_.topoChanger())
    )
{}


Foam::functionObjects::meshToMeshAdjustTimeStep::meshToMeshAdjustTimeStep
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    meshToMesh_
    (
        refCast<const fvMeshTopoChangers::meshToMesh>(mesh_.topoChanger())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::meshToMeshAdjustTimeStep::~meshToMeshAdjustTimeStep()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::meshToMeshAdjustTimeStep::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


Foam::scalar Foam::functionObjects::meshToMeshAdjustTimeStep::timeToNextAction()
{
    return meshToMesh_.timeToNextMesh();
}


bool Foam::functionObjects::meshToMeshAdjustTimeStep::execute()
{
    return true;
}


bool Foam::functionObjects::meshToMeshAdjustTimeStep::write()
{
    return true;
}


// ************************************************************************* //
