/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "faceZoneSelection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace faceSelections
{
    defineTypeNameAndDebug(faceZoneSelection, 0);
    addToRunTimeSelectionTable(faceSelection, faceZoneSelection, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceSelections::faceZoneSelection::faceZoneSelection
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    faceSelection(name, mesh, dict),
    zoneName_(dict_.lookup("zoneName"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceSelections::faceZoneSelection::~faceZoneSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceSelections::faceZoneSelection::select
(
    const label zoneID,
    labelList& faceToZoneID,
    boolList& faceToFlip
) const
{
    label readID = mesh_.faceZones().findZoneID(zoneName_);

    if (readID == -1)
    {
        FatalErrorIn
        (
            "faceSelections::faceZoneSelection::select(labelList&) const"
        )   << "Cannot find faceZone " << zoneName_ << nl << "Valid zones are "
            << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[readID];

    forAll(fZone, i)
    {
        label faceI = fZone[i];

        if (faceToZoneID[faceI] == -1)
        {
            faceToZoneID[faceI] = zoneID;
            faceToFlip[faceI] = fZone.flipMap()[i];
        }
        else if (faceToZoneID[faceI] != zoneID)
        {
            FatalErrorIn
            (
                "faceSelections::faceZoneSelection::select(labelList&) const"
            )   << "Face " << faceI << " already in faceZone "
                << faceToZoneID[faceI]
                << exit(FatalError);
        }
    }

    faceSelection::select(zoneID, faceToZoneID, faceToFlip);
}


// ************************************************************************* //
