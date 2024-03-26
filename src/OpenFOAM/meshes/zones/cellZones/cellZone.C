/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "cellZone.H"
#include "cellZones.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef Zone<cellZone, cellZones> cellZoneType;
    defineTemplateRunTimeSelectionTable(cellZoneType, dictionary);

    defineTypeNameAndDebug(cellZone, 0);
    addToRunTimeSelectionTable(cellZone, cellZone, dictionary);
}

const char * const Foam::cellZone::labelsName = "cellLabels";


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellZone::whichCell(const label globalCellID) const
{
    return Zone<cellZone, cellZones>::localIndex(globalCellID);
}


bool Foam::cellZone::checkDefinition(const bool report) const
{
    return Zone<cellZone, cellZones>::checkDefinition
    (
        zones_.mesh().nCells(),
        report
    );
}


void Foam::cellZone::topoChange(const polyTopoChangeMap& map)
{
    clearAddressing();

    labelHashSet indices;
    const labelList& cellMap = map.cellMap();
    const labelList& reverseCellMap = map.reverseCellMap();

    forAll(cellMap, celli)
    {
        if (cellMap[celli] >= 0 && localIndex(cellMap[celli]) != -1)
        {
            indices.insert(celli);
        }
    }

    forAll(reverseCellMap, celli)
    {
        if (reverseCellMap[celli] >= 0 && localIndex(celli) != -1)
        {
            indices.insert(reverseCellMap[celli]);
        }
    }

    labelList::operator=(indices.sortedToc());
}


void Foam::cellZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry(os, this->labelsName, *this);

    os  << token::END_BLOCK << endl;
}


// ************************************************************************* //
