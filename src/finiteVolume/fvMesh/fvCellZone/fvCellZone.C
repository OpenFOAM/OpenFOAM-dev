/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "fvCellZone.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvCellZone::setV()
{
    const labelList& cells(this->zone());

    V_ = 0;
    forAll(cells, i)
    {
        V_ += mesh_.V()[cells[i]];
    }
    reduce(V_, sumOp<scalar>());
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvCellZone::writeFileHeader
(
    const functionObjects::writeFile& wf,
    Ostream& file
)
{
    wf.writeCommented(file, "Selection");
    file<< setw(1) << ':' << setw(1) << ' ' << zoneName() << endl;
    wf.writeHeaderValue(file, "Volume", V());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvCellZone::fvCellZone(const fvMesh& mesh)
:
    generatedCellZone(mesh),
    mesh_(mesh),
    V_(gSum(mesh_.V()))
{}


Foam::fvCellZone::fvCellZone(const fvMesh& mesh, const dictionary& dict)
:
    generatedCellZone(mesh, dict),
    mesh_(mesh),
    V_(NaN)
{
    setV();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvCellZone::~fvCellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvCellZone::movePoints()
{
    generatedCellZone::movePoints();
    setV();
}


void Foam::fvCellZone::topoChange(const polyTopoChangeMap& map)
{
    generatedCellZone::topoChange(map);
    setV();
}


void Foam::fvCellZone::mapMesh(const polyMeshMap& map)
{
    generatedCellZone::mapMesh(map);
    setV();
}


void Foam::fvCellZone::distribute(const polyDistributionMap& map)
{
    generatedCellZone::distribute(map);
    setV();
}


bool Foam::fvCellZone::read(const dictionary& dict)
{
    generatedCellZone::read(dict);
    setV();

    return true;
}


// ************************************************************************* //
