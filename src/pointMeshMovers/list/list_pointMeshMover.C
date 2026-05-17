/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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

#include "list_pointMeshMover.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(list, 0);
    addToRunTimeSelectionTable
    (
        pointMeshMover,
        list,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::list::list
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointMeshMover(mesh, typeName),
    movers_(0)
{
    const dictionary& solversDict = dict.optionalSubDict("list");

    forAllConstIter(dictionary, solversDict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict = iter().dict();

            movers_.append
            (
                name,
                pointMeshMover::New(mesh, dict).ptr()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::list::~list()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::pointMeshMovers::list::newPoints()
{
    if (movers_.size())
    {
        // Accumulated displacement
        pointField disp(poly().nPoints(), Zero);

        forAllIter(PtrListDictionary<pointMeshMover>, movers_, iter)
        {
            disp += iter().newPoints() - poly().points();
        }

        return poly().points() + disp;
    }
    else
    {
        return poly().points();
    }
}


void Foam::pointMeshMovers::list::topoChange(const polyTopoChangeMap& map)
{
    forAllIter(PtrListDictionary<pointMeshMover>, movers_, iter)
    {
        iter().topoChange(map);
    }
}


void Foam::pointMeshMovers::list::mapMesh(const polyMeshMap& map)
{
    forAllIter(PtrListDictionary<pointMeshMover>, movers_, iter)
    {
        iter().mapMesh(map);
    }
}


void Foam::pointMeshMovers::list::distribute
(
    const polyDistributionMap& map
)
{
    forAllIter(PtrListDictionary<pointMeshMover>, movers_, iter)
    {
        iter().distribute(map);
    }
}


void Foam::pointMeshMovers::list::movePoints(const pointField& points)
{
    forAllIter(PtrListDictionary<pointMeshMover>, movers_, iter)
    {
        iter().movePoints(points);
    }
}


// ************************************************************************* //
