/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

#include "pointMeshMover_fvMeshMover.H"
#include "pointMeshMover.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(pointMeshMover, 0);
    addToRunTimeSelectionTable(fvMeshMover, pointMeshMover, fvMesh);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvMeshMover,
        pointMeshMover,
        fvMesh,
        motionSolver,
        "motionSolver"
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::pointMeshMover::pointMeshMover
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshMover(mesh),
    pointMeshMoverName_
    (
        dict.found("motionSolver") ? "motionSolver" : pointMeshMover::typeName
    ),
    pointMeshMoverPtr_
    (
        Foam::pointMeshMover::New
        (
            mesh,
            dict.optionalTypeDict(pointMeshMoverName_)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::pointMeshMover::~pointMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointMeshMover& Foam::fvMeshMovers::pointMeshMover::mover() const
{
    return pointMeshMoverPtr_();
}


bool Foam::fvMeshMovers::pointMeshMover::solidBodyMotion() const
{
    return pointMeshMoverPtr_->solidBodyMotion();
}


bool Foam::fvMeshMovers::pointMeshMover::update()
{
    mesh().preChange();

    mesh().movePoints(pointMeshMoverPtr_->newPoints());

    return true;
}


void Foam::fvMeshMovers::pointMeshMover::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMoverPtr_->topoChange(map);
}


void Foam::fvMeshMovers::pointMeshMover::mapMesh(const polyMeshMap& map)
{
    pointMeshMoverPtr_->mapMesh(map);
}


void Foam::fvMeshMovers::pointMeshMover::distribute
(
    const polyDistributionMap& map
)
{
    pointMeshMoverPtr_->distribute(map);
}


bool Foam::fvMeshMovers::pointMeshMover::write(const bool write) const
{
    if (write)
    {
        return pointMeshMoverPtr_->write();
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
