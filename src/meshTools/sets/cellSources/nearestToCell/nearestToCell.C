/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "nearestToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(nearestToCell, 0);

addToRunTimeSelectionTable(topoSetSource, nearestToCell, word);

addToRunTimeSelectionTable(topoSetSource, nearestToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::nearestToCell::usage_
(
    nearestToCell::typeName,
    "\n    Usage: nearestToCell (pt0 .. ptn)\n\n"
    "    Select the nearest cell for each of the points pt0 ..ptn\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearestToCell::combine(topoSet& set, const bool add) const
{
    forAll(points_, pointI)
    {
        addOrDelete(set, mesh_.findNearestCell(points_[pointI]), add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    const pointField& points
)
:
    topoSetSource(mesh),
    points_(points)
{}


// Construct from dictionary
Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    points_(dict.lookup("points"))
{}


// Construct from Istream
Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    points_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearestToCell::~nearestToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearestToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells nearest to " << points_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells nearest to " << points_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
