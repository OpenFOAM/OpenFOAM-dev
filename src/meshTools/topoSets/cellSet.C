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

#include "cellSet.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellSet, 0);

    addToRunTimeSelectionTable(topoSet, cellSet, word);
    addToRunTimeSelectionTable(topoSet, cellSet, size);
    addToRunTimeSelectionTable(topoSet, cellSet, set);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::cellSet::cellSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    // Make sure set within valid range
    check(mesh.nCells());
}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSet::~cellSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void Foam::cellSet::topoChange(const polyTopoChangeMap& map)
{
    updateLabels(map.reverseCellMap());
}


void Foam::cellSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.cellCentres(), maxLen);
}


// ************************************************************************* //
