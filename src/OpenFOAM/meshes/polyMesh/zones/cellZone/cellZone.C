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
#include "meshCellZones.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellZone, 0);
    defineRunTimeSelectionTable(cellZone, dictionary);
    addToRunTimeSelectionTable(cellZone, cellZone, dictionary);
}

const char * const Foam::cellZone::labelsName = "cellLabels";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellZone::cellZone
(
    const word& name,
    const labelUList& addr,
    const meshCellZones& mz
)
:
    zone(name, addr),
    meshZones_(mz)
{}


Foam::cellZone::cellZone
(
    const word& name,
    labelList&& addr,
    const meshCellZones& mz
)
:
    zone(name, move(addr)),
    meshZones_(mz)
{}


Foam::cellZone::cellZone
(
    const word& name,
    const dictionary& dict,
    const meshCellZones& mz
)
:
    zone(name, dict, this->labelsName),
    meshZones_(mz)
{}


Foam::cellZone::cellZone
(
    const cellZone& cz,
    const labelUList& addr,
    const meshCellZones& mz
)
:
    zone(cz, addr),
    meshZones_(mz)
{}


Foam::cellZone::cellZone
(
    const cellZone& cz,
    labelList&& addr,
    const meshCellZones& mz
)
:
    zone(cz, move(addr)),
    meshZones_(mz)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellZone::~cellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellZone::whichCell(const label globalCellID) const
{
    return zone::localIndex(globalCellID);
}


const Foam::meshCellZones& Foam::cellZone::meshZones() const
{
    return meshZones_;
}


bool Foam::cellZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(meshZones_.mesh().nCells(), report);
}


void Foam::cellZone::topoChange(const polyTopoChangeMap& map)
{
    clearAddressing();

    labelHashSet newIndices;
    const labelList& cellMap = map.cellMap();

    forAll(cellMap, celli)
    {
        if (cellMap[celli] >= 0 && localIndex(cellMap[celli]) != -1)
        {
            newIndices.insert(celli);
        }
    }

    labelList::operator=(newIndices.sortedToc());
}


void Foam::cellZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry(os, this->labelsName, *this);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cellZone::operator=(const cellZone& zn)
{
    zone::operator=(zn);
}


void Foam::cellZone::operator=(cellZone&& zn)
{
    zone::operator=(move(zn));
}


// ************************************************************************* //
