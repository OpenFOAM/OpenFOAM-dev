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

void Foam::fvCellZone::update()
{
    const labelList& cells(this->zone());

    nGlobalCells_ = cells.size();
    reduce(nGlobalCells_, sumOp<label>());

    V_ = 0;
    forAll(cells, i)
    {
        V_ += mesh_.V()[cells[i]];
    }
    reduce(V_, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvCellZone::fvCellZone(const fvMesh& mesh)
:
    generatedCellZone(mesh),
    mesh_(mesh),
    nGlobalCells_(returnReduce(mesh.nCells(), sumOp<label>())),
    V_(gSum(mesh_.V()))
{}


Foam::fvCellZone::fvCellZone(const fvMesh& mesh, const dictionary& dict)
:
    generatedCellZone(mesh, dict),
    mesh_(mesh),
    nGlobalCells_(-1),
    V_(NaN)
{
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvCellZone::~fvCellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvCellZone::writeFileHeader
(
    const functionObjects::writeFile& wf,
    Ostream& file
)
{
    wf.writeCommented(file, "Selection");
    file<< setw(1) << ':' << setw(1) << ' ' << name() << endl;
    wf.writeHeaderValue(file, "Cells", nGlobalCells());
    wf.writeHeaderValue(file, "Volume", V());
}


void Foam::fvCellZone::movePoints()
{
    generatedCellZone::movePoints();
    update();
}


void Foam::fvCellZone::topoChange(const polyTopoChangeMap& map)
{
    generatedCellZone::topoChange(map);
    update();
}


void Foam::fvCellZone::mapMesh(const polyMeshMap& map)
{
    generatedCellZone::mapMesh(map);
    update();
}


void Foam::fvCellZone::distribute(const polyDistributionMap& map)
{
    generatedCellZone::distribute(map);
    update();
}


bool Foam::fvCellZone::read(const dictionary& dict)
{
    generatedCellZone::read(dict);
    update();

    return true;
}


// ************************************************************************* //
