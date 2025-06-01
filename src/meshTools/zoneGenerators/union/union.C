/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "union.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(Union, 0);
        addToRunTimeSelectionTable(zoneGenerator, Union, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneGenerators::Union::select
(
    boolList& selected,
    const labelList& indices,
    const bool select
) const
{
    forAll(indices, i)
    {
        selected[indices[i]] = select;
    }
}


void Foam::zoneGenerators::Union::select
(
    boolList& selected,
    boolList& flipMap,
    const labelList& indices,
    const boolList& zoneFlipMap,
    const bool select
) const
{
    forAll(indices, i)
    {
        selected[indices[i]] = select;
        flipMap[indices[i]] = zoneFlipMap[i];
    }
}


Foam::zoneSet Foam::zoneGenerators::Union::generate
(
    const bool diff,
    const bool all
) const
{
    boolList selectedPoints;
    boolList selectedCells;
    boolList selectedFaces;
    boolList flipMap;

    forAll(zoneGenerators_, i)
    {
        zoneSet zs(zoneGenerators_[i].generate());

        // Select or deselect
        const bool sel = all ? false : (diff ? i == 0 : true);

        if
        (
            zoneType_ == zoneTypesAll::point
         || (zoneType_ == zoneTypesAll::all && zs.pZone.valid())
        )
        {
            selectedPoints.setSize(mesh_.nPoints(), all);
            select(selectedPoints, zs.pZone, sel);
        }

        if
        (
            zoneType_ == zoneTypesAll::cell
         || (zoneType_ == zoneTypesAll::all && zs.cZone.valid())
        )
        {
            selectedCells.setSize(mesh_.nCells(), all);
            select(selectedCells, zs.cZone, sel);
        }

        if
        (
            zoneType_ == zoneTypesAll::face
         || (zoneType_ == zoneTypesAll::all && zs.fZone.valid())
        )
        {
            selectedFaces.setSize(mesh_.nFaces(), all);
            flipMap.setSize(mesh_.nFaces(), false);
            select
            (
                selectedFaces,
                flipMap,
                zs.fZone(),
                zs.fZone().flipMap(),
                sel
            );
        }
    }

    moveUpdate_ = zoneGenerators_.moveUpdate();

    const labelList faceIndices(indices(selectedFaces));

    return zoneSet
    (
        selectedPoints.size()
      ? new pointZone
        (
            zoneName_,
            indices(selectedPoints),
            mesh_.pointZones(),
            moveUpdate_,
            true
        )
      : nullptr,
        selectedCells.size()
      ? new cellZone
        (
            zoneName_,
            indices(selectedCells),
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
      : nullptr,
        selectedFaces.size()
      ? new faceZone
        (
            zoneName_,
            faceIndices,
            boolList(flipMap, faceIndices),
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
      : nullptr
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::Union::Union
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneType_
    (
        zoneTypesAllNames.lookupOrDefault("zoneType", dict, zoneTypesAll::all)
    ),
    zoneGenerators_(mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::Union::~Union()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::Union::generate() const
{
    return generate(false, false);
}


// ************************************************************************* //
