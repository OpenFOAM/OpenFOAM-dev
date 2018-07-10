/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "faceToCell.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, faceToCell, word);
    addToRunTimeSelectionTable(topoSetSource, faceToCell, istream);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::faceToCell::faceAction,
        4
    >::names[] =
    {
        "neighbour",
        "owner",
        "any",
        "all"
    };
}


Foam::topoSetSource::addToUsageTable Foam::faceToCell::usage_
(
    faceToCell::typeName,
    "\n    Usage: faceToCell <faceSet> neighbour|owner|any|all\n\n"
    "    Select cells that are the owner|neighbour|any"
    " of the faces in the faceSet or where all faces are in the faceSet\n\n"
);

const Foam::NamedEnum<Foam::faceToCell::faceAction, 4>
    Foam::faceToCell::faceActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceToCell::combine(topoSet& set, const bool add) const
{
    // Load the set
    faceSet loadedSet(mesh_, setName_);


    // Handle owner/neighbour/any selection
    forAllConstIter(faceSet, loadedSet, iter)
    {
        const label facei = iter.key();

        if ((option_ == OWNER) || (option_ == ANY))
        {
            const label celli = mesh_.faceOwner()[facei];

            addOrDelete(set, celli, add);
        }

        if (mesh_.isInternalFace(facei))
        {
            if ((option_ == NEIGHBOUR) || (option_ == ANY))
            {
                const label celli = mesh_.faceNeighbour()[facei];

                addOrDelete(set, celli, add);
            }
        }
    }

    // Handle all selection.
    if (option_ == ALL)
    {
        // Count number of selected faces per cell.

        Map<label> facesPerCell(loadedSet.size());

        forAllConstIter(faceSet, loadedSet, iter)
        {
            const label facei = iter.key();
            const label own = mesh_.faceOwner()[facei];

            Map<label>::iterator fndOwn = facesPerCell.find(own);

            if (fndOwn == facesPerCell.end())
            {
                facesPerCell.insert(own, 1);
            }
            else
            {
                fndOwn()++;
            }

            if (mesh_.isInternalFace(facei))
            {
                label nei = mesh_.faceNeighbour()[facei];

                Map<label>::iterator fndNei = facesPerCell.find(nei);

                if (fndNei == facesPerCell.end())
                {
                    facesPerCell.insert(nei, 1);
                }
                else
                {
                    fndNei()++;
                }
            }
        }

        // Include cells that are referenced as many times as they have faces
        // -> all faces in set.
        forAllConstIter(Map<label>, facesPerCell, iter)
        {
            const label celli = iter.key();

            if (iter() == mesh_.cells()[celli].size())
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceToCell::faceToCell
(
    const polyMesh& mesh,
    const word& setName,
    const faceAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


Foam::faceToCell::faceToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(faceActionNames_.read(dict.lookup("option")))
{}


Foam::faceToCell::faceToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceToCell::~faceToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells according to faceSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells according to faceSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
