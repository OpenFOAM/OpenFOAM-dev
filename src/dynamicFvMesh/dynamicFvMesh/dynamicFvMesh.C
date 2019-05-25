/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "dynamicFvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFvMesh, 0);
    defineRunTimeSelectionTable(dynamicFvMesh, IOobject);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::dynamicFvMesh::dynamicMeshDictIOobject(const IOobject& io)
{
    // defaultRegion (region0) gets loaded from constant, other ones get loaded
    // from constant/<regionname>. Normally we'd use polyMesh::dbDir() but we
    // haven't got a polyMesh yet ...
    return IOobject
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::dynamicFvMesh(const IOobject& io)
:
    fvMesh(io),
    dynamicMeshDict_(IOdictionary(dynamicMeshDictIOobject(io)))
{}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        move(points),
        move(faces),
        move(allOwner),
        move(allNeighbour),
        syncPar
    ),
    dynamicMeshDict_(IOdictionary(dynamicMeshDictIOobject(io)))
{}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        move(points),
        move(faces),
        move(cells),
        syncPar
    ),
    dynamicMeshDict_(IOdictionary(dynamicMeshDictIOobject(io)))
{}


Foam::dynamicFvMesh::velocityMotionCorrection::velocityMotionCorrection
(
    const dynamicFvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    velocityFields_(dict.lookupOrDefault("velocityFields", wordList()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::~dynamicFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicFvMesh::velocityMotionCorrection::update() const
{
    forAll(velocityFields_, i)
    {
        if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
        {
            mesh_.lookupObjectRef<volVectorField>
            (
                velocityFields_[i]
            ).correctBoundaryConditions();
        }
    }
}

// ************************************************************************* //
