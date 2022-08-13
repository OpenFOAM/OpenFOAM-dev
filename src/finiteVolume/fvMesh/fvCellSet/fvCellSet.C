/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fvCellSet.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvCellSet::setV()
{
    Info<< incrIndent;

    const labelList& cells = this->cells();

    V_ = 0;
    forAll(cells, i)
    {
        V_ += mesh_.V()[cells[i]];
    }
    reduce(V_, sumOp<scalar>());

    Info<< indent
        << "- selected " << returnReduce(cells.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << endl;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvCellSet::fvCellSet
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    polyCellSet(mesh, dict),
    mesh_(mesh),
    V_(NaN)
{
    setV();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvCellSet::~fvCellSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvCellSet::movePoints()
{
    polyCellSet::movePoints();
    setV();
}


void Foam::fvCellSet::topoChange(const polyTopoChangeMap& map)
{
    polyCellSet::topoChange(map);
    setV();
}


void Foam::fvCellSet::mapMesh(const polyMeshMap& map)
{
    polyCellSet::mapMesh(map);
    setV();
}


void Foam::fvCellSet::distribute(const polyDistributionMap& map)
{
    polyCellSet::distribute(map);
    setV();
}


bool Foam::fvCellSet::read(const dictionary& dict)
{
    polyCellSet::read(dict);
    setV();

    return true;
}


// ************************************************************************* //
