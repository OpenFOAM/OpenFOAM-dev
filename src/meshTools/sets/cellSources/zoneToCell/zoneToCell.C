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

#include "zoneToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToCell::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    forAll(mesh_.cellZones(), i)
    {
        const cellZone& zone = mesh_.cellZones()[i];

        if (zoneName_.match(zone.name()))
        {
            const labelList& cellLabels = mesh_.cellZones()[i];

            Info<< "    Found matching zone " << zone.name()
                << " with " << returnReduce(cellLabels.size(), sumOp<label>())
                << " cells." << endl;

            hasMatched = true;

            forAll(cellLabels, i)
            {
                // Only do active cells
                if (cellLabels[i] < mesh_.nCells())
                {
                    addOrDelete(set, cellLabels[i], add);
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningInFunction
            << "Cannot find any cellZone named " << zoneName_ << endl
            << "Valid names are " << mesh_.cellZones().names() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const word& zoneName
)
:
    topoSetSource(mesh),
    zoneName_(zoneName)
{}


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    zoneName_
    (
        dict.dictName() == "sourceInfo"
      ? dict.lookupBackwardsCompatible({"zone", "name"})
      : dict.lookup("zone")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneToCell::~zoneToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells of cellZone " << zoneName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of cellZone " << zoneName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
