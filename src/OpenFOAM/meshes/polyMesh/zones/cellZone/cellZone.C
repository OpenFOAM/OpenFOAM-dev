/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "cellZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "IOstream.H"
#include "demandDrivenData.H"

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
    const label index,
    const cellZoneMesh& zm
)
:
    zone(name, addr, index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const word& name,
    labelList&& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(name, move(addr), index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(name, dict, this->labelsName, index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const cellZone& cz,
    const labelUList& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(cz, addr, index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const cellZone& cz,
    labelList&& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(cz, move(addr), index),
    zoneMesh_(zm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellZone::~cellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellZone::whichCell(const label globalCellID) const
{
    return zone::localID(globalCellID);
}


const Foam::cellZoneMesh& Foam::cellZone::zoneMesh() const
{
    return zoneMesh_;
}


bool Foam::cellZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(zoneMesh_.mesh().nCells(), report);
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
    clearAddressing();
    zone::operator=(zn);
}


void Foam::cellZone::operator=(cellZone&& zn)
{
    clearAddressing();
    zone::operator=(move(zn));
}


void Foam::cellZone::operator=(const labelUList& addr)
{
    clearAddressing();
    zone::operator=(addr);
}


void Foam::cellZone::operator=(labelList&& addr)
{
    clearAddressing();
    zone::operator=(move(addr));
}


// ************************************************************************* //
