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

#include "setToCellZone.H"
#include "polyMesh.H"
#include "cellZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setToCellZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setToCellZone, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setToCellZone::setToCellZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetSource(mesh),
    setName_(setName)
{}


Foam::setToCellZone::setToCellZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookupBackwardsCompatible<word>({"set", "cellSet"}))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setToCellZone::~setToCellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToCellZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<cellZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a cellZoneSet." << endl;
    }
    else
    {
        cellZoneSet& fzSet = refCast<cellZoneSet>(set);

        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding all cells from cellSet " << setName_
                << " ..." << endl;

            // Load the sets
            cellSet fSet(mesh_, setName_);

            // Start off from copy
            DynamicList<label> newAddressing(fzSet.addressing());

            forAllConstIter(cellSet, fSet, iter)
            {
                label celli = iter.key();

                if (!fzSet.found(celli))
                {
                    newAddressing.append(celli);
                }
            }

            fzSet.addressing().transfer(newAddressing);
            fzSet.updateSet();
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing all cells from cellSet " << setName_
                << " ..." << endl;

            // Load the set
            cellSet loadedSet(mesh_, setName_);

            // Start off empty
            DynamicList<label> newAddressing(fzSet.addressing().size());

            forAll(fzSet.addressing(), i)
            {
                if (!loadedSet.found(fzSet.addressing()[i]))
                {
                    newAddressing.append(fzSet.addressing()[i]);
                }
            }
            fzSet.addressing().transfer(newAddressing);
            fzSet.updateSet();
        }
    }
}


// ************************************************************************* //
