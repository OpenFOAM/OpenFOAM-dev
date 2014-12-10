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

#include "zoneToFace.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(zoneToFace, 0);

addToRunTimeSelectionTable(topoSetSource, zoneToFace, word);

addToRunTimeSelectionTable(topoSetSource, zoneToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::zoneToFace::usage_
(
    zoneToFace::typeName,
    "\n    Usage: zoneToFace zone\n\n"
    "    Select all faces in the faceZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToFace::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    forAll(mesh_.faceZones(), i)
    {
        const faceZone& zone = mesh_.faceZones()[i];

        if (zoneName_.match(zone.name()))
        {
            const labelList& faceLabels = mesh_.faceZones()[i];

            Info<< "    Found matching zone " << zone.name()
                << " with " << faceLabels.size() << " faces." << endl;

            hasMatched = true;

            forAll(faceLabels, i)
            {
                // Only do active faces
                if (faceLabels[i] < mesh_.nFaces())
                {
                    addOrDelete(set, faceLabels[i], add);
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningIn("zoneToFace::combine(topoSet&, const bool)")
            << "Cannot find any faceZone named " << zoneName_ << endl
            << "Valid names are " << mesh_.faceZones().names() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const word& zoneName
)
:
    topoSetSource(mesh),
    zoneName_(zoneName)
{}


// Construct from dictionary
Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    zoneName_(dict.lookup("name"))
{}


// Construct from Istream
Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    zoneName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneToFace::~zoneToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all faces of faceZone " << zoneName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all faces of faceZone " << zoneName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
