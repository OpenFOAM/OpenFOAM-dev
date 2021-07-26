/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "setsToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "cellSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setsToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setsToFaceZone, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setsToFaceZone::setsToFaceZone
(
    const polyMesh& mesh,
    const word& faceSetName,
    const word& cellSetName,
    const Switch& flip
)
:
    topoSetSource(mesh),
    faceSetName_(faceSetName),
    cellSetName_(cellSetName),
    flip_(flip)
{}


Foam::setsToFaceZone::setsToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    faceSetName_(dict.lookup("faceSet")),
    cellSetName_(dict.lookup("cellSet")),
    flip_(dict.lookupOrDefault("flip", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setsToFaceZone::~setsToFaceZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setsToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet." << endl;
    }
    else
    {
        faceZoneSet& fzSet = refCast<faceZoneSet>(set);

        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding all faces from faceSet " << faceSetName_
                << " ..." << endl;

            // Load the sets
            faceSet fSet(mesh_, faceSetName_);
            cellSet cSet(mesh_, cellSetName_);

            // Start off from copy
            DynamicList<label> newAddressing(fzSet.addressing());
            DynamicList<bool> newFlipMap(fzSet.flipMap());

            forAllConstIter(faceSet, fSet, iter)
            {
                label facei = iter.key();

                if (!fzSet.found(facei))
                {
                    bool flipFace = false;

                    label own = mesh_.faceOwner()[facei];
                    bool ownFound = cSet.found(own);

                    if (mesh_.isInternalFace(facei))
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        bool neiFound = cSet.found(nei);

                        if (ownFound && !neiFound)
                        {
                            flipFace = false;
                        }
                        else if (!ownFound && neiFound)
                        {
                            flipFace = true;
                        }
                        else
                        {
                            WarningInFunction
                                << "One of owner or neighbour of internal face "
                                << facei << " should be in cellSet "
                                << cSet.name()
                                << " to be able to determine orientation."
                                << endl
                                << "Face:" << facei << " own:" << own
                                << " OwnInCellSet:" << ownFound
                                << " nei:" << nei
                                << " NeiInCellSet:" << neiFound
                                << endl;
                        }
                    }
                    else
                    {
                        flipFace = !ownFound;
                    }


                    if (flip_)
                    {
                        flipFace = !flipFace;
                    }

                    newAddressing.append(facei);
                    newFlipMap.append(flipFace);
                }
            }

            fzSet.addressing().transfer(newAddressing);
            fzSet.flipMap().transfer(newFlipMap);
            fzSet.updateSet();
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing all faces from faceSet " << faceSetName_
                << " ..." << endl;

            // Load the set
            faceZoneSet loadedSet(mesh_, faceSetName_);

            // Start off empty
            DynamicList<label> newAddressing(fzSet.addressing().size());
            DynamicList<bool> newFlipMap(fzSet.flipMap().size());

            forAll(fzSet.addressing(), i)
            {
                if (!loadedSet.found(fzSet.addressing()[i]))
                {
                    newAddressing.append(fzSet.addressing()[i]);
                    newFlipMap.append(fzSet.flipMap()[i]);
                }
            }
            fzSet.addressing().transfer(newAddressing);
            fzSet.flipMap().transfer(newFlipMap);
            fzSet.updateSet();
        }
    }
}


// ************************************************************************* //
