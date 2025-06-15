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

#include "pointSet.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointSet, 0);

    addToRunTimeSelectionTable(topoSet, pointSet, word);
    addToRunTimeSelectionTable(topoSet, pointSet, size);
    addToRunTimeSelectionTable(topoSet, pointSet, set);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointSet::pointSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    check(mesh.nPoints());
}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::pointSet::pointSet
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

Foam::pointSet::~pointSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointSet::sync(const polyMesh& mesh)
{
    // Convert to boolList

    boolList contents(mesh.nPoints(), false);

    forAllConstIter(pointSet, *this, iter)
    {
        contents[iter.key()] = true;
    }
    syncTools::syncPointList
    (
        mesh,
        contents,
        orEqOp<bool>(),
        false           // null value
    );

    // Convert back to labelHashSet

    labelHashSet newContents(size());

    forAll(contents, pointi)
    {
        if (contents[pointi])
        {
            newContents.insert(pointi);
        }
    }

    transfer(newContents);
}


Foam::label Foam::pointSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nPoints();
}


void Foam::pointSet::topoChange(const polyTopoChangeMap& map)
{
    updateLabels(map.reversePointMap());
}


void Foam::pointSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.points(), maxLen);
}


// ************************************************************************* //
