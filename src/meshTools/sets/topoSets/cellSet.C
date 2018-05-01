/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cellSet, 0);

addToRunTimeSelectionTable(topoSet, cellSet, word);
addToRunTimeSelectionTable(topoSet, cellSet, size);
addToRunTimeSelectionTable(topoSet, cellSet, set);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cellSet::cellSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


cellSet::cellSet
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


cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


// Database constructors (for when no mesh available)
cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet
    (
        IOobject
        (
            name,
            runTime.findInstance
            (
                polyMesh::meshSubDir/"sets",    // polyMesh::meshSubDir,
                word::null,                     //"faces"
                IOobject::MUST_READ,
                runTime.findInstance
                (
                    polyMesh::meshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                )
            ),
            polyMesh::meshSubDir/"sets",
            runTime,
            r,
            w
        ),
        typeName
    )
{}


cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet
    (
        IOobject
        (
            name,
            runTime.findInstance
            (
                polyMesh::meshSubDir/"sets",    // polyMesh::meshSubDir,
                word::null,                     //"faces"
                IOobject::NO_READ,
                runTime.findInstance
                (
                    polyMesh::meshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                )
            ),
            polyMesh::meshSubDir/"sets",
            runTime,
            IOobject::NO_READ,
            w
        ),
        size
    )
{}


cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    topoSet
    (
        IOobject
        (
            name,
            runTime.findInstance
            (
                polyMesh::meshSubDir/"sets",    // polyMesh::meshSubDir,
                word::null,                     //"faces"
                IOobject::NO_READ,
                runTime.findInstance
                (
                    polyMesh::meshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                )
            ),
            polyMesh::meshSubDir/"sets",
            runTime,
            IOobject::NO_READ,
            w
        ),
        set
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellSet::~cellSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label cellSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void cellSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseCellMap());
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
