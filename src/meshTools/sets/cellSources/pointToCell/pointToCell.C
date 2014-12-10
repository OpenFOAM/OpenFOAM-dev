/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "pointToCell.H"
#include "polyMesh.H"
#include "pointSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, pointToCell, word);
    addToRunTimeSelectionTable(topoSetSource, pointToCell, istream);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::pointToCell::pointAction,
        2
    >::names[] =
    {
        "any",
        "edge"
    };
}


Foam::topoSetSource::addToUsageTable Foam::pointToCell::usage_
(
    pointToCell::typeName,
    "\n    Usage: pointToCell <pointSet> any|edge\n\n"
    "    Select all cells with any point ('any') or any edge ('edge')"
    " in the pointSet\n\n"
);

const Foam::NamedEnum<Foam::pointToCell::pointAction, 2>
    Foam::pointToCell::pointActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointToCell::combine(topoSet& set, const bool add) const
{
    // Load the set
    pointSet loadedSet(mesh_, setName_);


    // Handle any selection
    if (option_ == ANY)
    {
        forAllConstIter(pointSet, loadedSet, iter)
        {
            const label pointI = iter.key();
            const labelList& pCells = mesh_.pointCells()[pointI];

            forAll(pCells, pCellI)
            {
                addOrDelete(set, pCells[pCellI], add);
            }
        }
    }
    else if (option_ == EDGE)
    {
        const faceList& faces = mesh_.faces();
        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            forAll(f, fp)
            {
                if (loadedSet.found(f[fp]) && loadedSet.found(f.nextLabel(fp)))
                {
                    addOrDelete(set, mesh_.faceOwner()[faceI], add);
                    if (mesh_.isInternalFace(faceI))
                    {
                        addOrDelete(set, mesh_.faceNeighbour()[faceI], add);
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


// Construct from dictionary
Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(pointActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointToCell::~pointToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
