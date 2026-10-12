/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "solidDisplacement_fvMeshMover.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(solidDisplacement, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        solidDisplacement,
        fvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::solidDisplacement::solidDisplacement
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshMover(mesh),
    pointMeshMovers::displacement(mesh.poly(), dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::solidDisplacement::~solidDisplacement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshMovers::solidDisplacement::update()
{
    mesh().preChange();

    mesh().movePoints(newPoints());

    return true;
}


Foam::tmp<Foam::pointField>
Foam::fvMeshMovers::solidDisplacement::newPoints()
{
    // The points have moved so before interpolation update
    // the pointMeshMover accordingly
    movePoints(mesh().points());

    pointDisplacement_.boundaryRef().updateCoeffs();

    // Interpolate the displacement volVectorField
    volPointInterpolation::New(mesh()).interpolate
    (
        mesh().lookupObject<volVectorField>("D"),
        pointDisplacement_
    );

    // Calculate the current points
    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );
    pointField& curPoints = tcurPoints.ref();

    // Constrain the current points for 2D cases
    twoDCorrectPoints(curPoints);

    return tcurPoints;
}


void Foam::fvMeshMovers::solidDisplacement::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::displacement::topoChange(map);
}


void Foam::fvMeshMovers::solidDisplacement::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::displacement::mapMesh(map);
}


// ************************************************************************* //
