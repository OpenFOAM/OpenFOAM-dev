/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "cellToFace.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "Time.H"
#include "syncTools.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToFace, word);
    addToRunTimeSelectionTable(topoSetSource, cellToFace, istream);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::cellToFace::cellAction,
        2
    >::names[] =
    {
        "all",
        "both"
    };
}


Foam::topoSetSource::addToUsageTable Foam::cellToFace::usage_
(
    cellToFace::typeName,
    "\n    Usage: cellToFace <cellSet> all|both\n\n"
    "    Select -all : all faces of cells in the cellSet\n"
    "           -both: faces where both neighbours are in the cellSet\n\n"
);

const Foam::NamedEnum<Foam::cellToFace::cellAction, 2>
    Foam::cellToFace::cellActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToFace::combine(topoSet& set, const bool add) const
{
    // Load the set
    if (!exists(mesh_.time().path()/topoSet::localPath(mesh_, setName_)))
    {
        SeriousError<< "Cannot load set "
            << setName_ << endl;
    }

    cellSet loadedSet(mesh_, setName_);

    if (option_ == ALL)
    {
        // Add all faces from cell
        forAllConstIter(cellSet, loadedSet, iter)
        {
            const label celli = iter.key();
            const labelList& cFaces = mesh_.cells()[celli];

            forAll(cFaces, cFacei)
            {
                addOrDelete(set, cFaces[cFacei], add);
            }
        }
    }
    else if (option_ == BOTH)
    {
        // Add all faces whose both neighbours are in set.

        const label nInt = mesh_.nInternalFaces();
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();


        // Check all internal faces
        for (label facei = 0; facei < nInt; facei++)
        {
            if (loadedSet.found(own[facei]) && loadedSet.found(nei[facei]))
            {
                addOrDelete(set, facei, add);
            }
        }


        // Get coupled cell status
        boolList neiInSet(mesh_.nFaces()-nInt, false);

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    neiInSet[facei-nInt] = loadedSet.found(own[facei]);
                    facei++;
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiInSet);


        // Check all boundary faces
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    if (loadedSet.found(own[facei]) && neiInSet[facei-nInt])
                    {
                        addOrDelete(set, facei, add);
                    }
                    facei++;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToFace::cellToFace
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


Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(cellActionNames_.read(dict.lookup("option")))
{}


Foam::cellToFace::cellToFace
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

Foam::cellToFace::~cellToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding faces according to cellSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces according to cellSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
