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

#include "cellToPoint.H"
#include "polyMesh.H"
#include "cellSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, cellToPoint, istream);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::cellToPoint::cellAction,
        1
    >::names[] =
    {
        "all"
    };
}


Foam::topoSetSource::addToUsageTable Foam::cellToPoint::usage_
(
    cellToPoint::typeName,
    "\n    Usage: cellToPoint <cellSet> all\n\n"
    "    Select all points of cells in the cellSet\n\n"
);

const Foam::NamedEnum<Foam::cellToPoint::cellAction, 1>
    Foam::cellToPoint::cellActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToPoint::combine(topoSet& set, const bool add) const
{
    // Load the set
    cellSet loadedSet(mesh_, setName_);

    // Add all point from cells in loadedSet
    forAllConstIter(cellSet, loadedSet, iter)
    {
        const label cellI = iter.key();
        const labelList& cFaces = mesh_.cells()[cellI];

        forAll(cFaces, cFaceI)
        {
            const face& f = mesh_.faces()[cFaces[cFaceI]];

            forAll(f, fp)
            {
                addOrDelete(set, f[fp], add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellToPoint::cellToPoint
(
    const polyMesh& mesh,
    const word& setName,
    const cellAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


// Construct from dictionary
Foam::cellToPoint::cellToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(cellActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::cellToPoint::cellToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(cellActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellToPoint::~cellToPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding from " << setName_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing from " << setName_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
