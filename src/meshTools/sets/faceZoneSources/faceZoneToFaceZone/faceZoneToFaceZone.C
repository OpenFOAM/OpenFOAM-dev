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

#include "faceZoneToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, faceZoneToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, faceZoneToFaceZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::faceZoneToFaceZone::usage_
(
    faceZoneToFaceZone::typeName,
    "\n    Usage: faceZoneToFaceZone <faceZone>\n\n"
    "    Select all faces in the faceZone\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneToFaceZone::faceZoneToFaceZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetSource(mesh),
    setName_(setName)
{}


Foam::faceZoneToFaceZone::faceZoneToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("zone"))
{}


Foam::faceZoneToFaceZone::faceZoneToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceZoneToFaceZone::~faceZoneToFaceZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZoneToFaceZone::applyToSet
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
        faceZoneSet& fSet = refCast<faceZoneSet>(set);

        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding all faces from faceZone " << setName_ << " ..."
                << endl;

            // Load the set
            faceZoneSet loadedSet(mesh_, setName_);

            DynamicList<label> newAddressing(fSet.addressing());
            DynamicList<bool> newFlipMap(fSet.flipMap());

            forAll(loadedSet.addressing(), i)
            {
                if (!fSet.found(loadedSet.addressing()[i]))
                {
                    newAddressing.append(loadedSet.addressing()[i]);
                    newFlipMap.append(loadedSet.flipMap()[i]);
                }
            }
            fSet.addressing().transfer(newAddressing);
            fSet.flipMap().transfer(newFlipMap);
            fSet.updateSet();
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing all faces from faceZone " << setName_ << " ..."
                << endl;

            // Load the set
            faceZoneSet loadedSet(mesh_, setName_);

            DynamicList<label> newAddressing(fSet.addressing().size());
            DynamicList<bool> newFlipMap(fSet.flipMap().size());

            forAll(fSet.addressing(), i)
            {
                if (!loadedSet.found(fSet.addressing()[i]))
                {
                    newAddressing.append(fSet.addressing()[i]);
                    newFlipMap.append(fSet.flipMap()[i]);
                }
            }
            fSet.addressing().transfer(newAddressing);
            fSet.flipMap().transfer(newFlipMap);
            fSet.updateSet();
        }
    }
}


// ************************************************************************* //
