/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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

#include "externalDisplacement_pointMeshMover.H"
#include "localPointRegion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(externalDisplacement, 0);

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        externalDisplacement,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::externalDisplacement::externalDisplacement
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointMeshMovers::displacement(mesh, dict, typeName),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::externalDisplacement::~externalDisplacement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::externalDisplacementMeshMover&
Foam::pointMeshMovers::externalDisplacement::meshMover() const
{
    if (!meshMoverPtr_.valid())
    {
        const word moverType(dict_.lookup("meshMover"));

        meshMoverPtr_ = externalDisplacementMeshMover::New
        (
            moverType,
            dict_,
            localPointRegion::findDuplicateFacePairs(poly()),
            pointDisplacement_
        );
    }
    return meshMoverPtr_();
}


Foam::tmp<Foam::pointField>
Foam::pointMeshMovers::externalDisplacement::newPoints()
{
    // The points have moved so before calculation update
    // the mesh and pointMeshMover accordingly
    movePoints(poly().points());

    // Update any point motion bcs (e.g. timevarying)
    pointDisplacement().boundaryFieldRef().updateCoeffs();

    label nAllowableErrors = 0;
    labelList checkFaces(identityMap(poly().nFaces()));
    meshMover().move
    (
        dict_,
        nAllowableErrors,
        checkFaces
    );

    // This will have updated the mesh and implicitly the pointDisplacement
    pointDisplacement().correctBoundaryConditions();

    return points();
}


void Foam::pointMeshMovers::externalDisplacement::movePoints
(
    const pointField& p
)
{
    pointMeshMovers::displacement::movePoints(p);

    // Update meshMover for new geometry
    if (meshMoverPtr_.valid())
    {
        meshMover().movePoints(p);
    }
}


void Foam::pointMeshMovers::externalDisplacement::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::displacement::topoChange(map);

    // Update meshMover for new topology
    meshMoverPtr_.clear();
}


// ************************************************************************* //
