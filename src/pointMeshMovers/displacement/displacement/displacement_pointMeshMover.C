/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "displacement_pointMeshMover.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(displacement, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::pointMeshMovers::displacement::points() const
{
    return points0() + pointDisplacement_.primitiveField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::displacement::displacement
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    displacementPoints0(mesh, dict, type),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::displacement::~displacement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointMeshMovers::displacement::mapMesh(const polyMeshMap& map)
{
    FatalErrorInFunction
        << "Mesh-to-mesh mapping in not implemented for displacement solvers"
        << nl
        << "    velocity based motion solvers are preferable for cases in which"
           " the mesh is reset periodically avoiding accumulation of error."
        << exit(FatalError);
}


// ************************************************************************* //
