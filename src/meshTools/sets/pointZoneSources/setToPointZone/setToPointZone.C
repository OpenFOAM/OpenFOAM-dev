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

#include "setToPointZone.H"
#include "polyMesh.H"
#include "pointZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setToPointZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setToPointZone, word);
    addToRunTimeSelectionTable(topoSetSource, setToPointZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::setToPointZone::usage_
(
    setToPointZone::typeName,
    "\n    Usage: setToPointZone <pointSet>\n\n"
    "    Select all points in the pointSet.\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setToPointZone::setToPointZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetSource(mesh),
    setName_(setName)
{}


Foam::setToPointZone::setToPointZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookupBackwardsCompatible<word>({"set", "pointSet"}))
{}


Foam::setToPointZone::setToPointZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setToPointZone::~setToPointZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToPointZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<pointZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a pointZoneSet." << endl;
    }
    else
    {
        pointZoneSet& fzSet = refCast<pointZoneSet>(set);

        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding all points from pointSet " << setName_
                << " ..." << endl;

            // Load the sets
            pointSet fSet(mesh_, setName_);

            // Start off from copy
            DynamicList<label> newAddressing(fzSet.addressing());

            forAllConstIter(pointSet, fSet, iter)
            {
                label pointi = iter.key();

                if (!fzSet.found(pointi))
                {
                    newAddressing.append(pointi);
                }
            }

            fzSet.addressing().transfer(newAddressing);
            fzSet.updateSet();
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing all points from pointSet " << setName_
                << " ..." << endl;

            // Load the set
            pointSet loadedSet(mesh_, setName_);

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
