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

#include "intersection.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(intersection, 0);
        addToRunTimeSelectionTable(zoneGenerator, intersection, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneGenerators::intersection::countIntersections
(
    labelList& nIntersections,
    const labelList& indices
) const
{
    forAll(indices, i)
    {
        nIntersections[indices[i]]++;
    }
}


void Foam::zoneGenerators::intersection::countIntersections
(
    labelList& nIntersections,
    boolList& flipMap,
    const labelList& indices,
    const boolList& zoneFlipMap
) const
{
    forAll(indices, i)
    {
        nIntersections[indices[i]]++;
        flipMap[indices[i]] = zoneFlipMap[i];
    }
}


Foam::labelList Foam::zoneGenerators::intersection::indices
(
    const labelList& nIntersections
) const
{
    const label nZones = zoneGenerators_.size();
    label nSelected = 0;
    forAll(nIntersections, i)
    {
        if (nIntersections[i] == nZones)
        {
            nSelected++;
        }
    }

    labelList intersectionIndices(nSelected);

    label ui = 0;
    forAll(nIntersections, i)
    {
        if (nIntersections[i] == nZones)
        {
            intersectionIndices[ui++] = i;
        }
    }

    return intersectionIndices;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::intersection::intersection
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
{
    moveUpdate_ = zoneGenerators_.moveUpdate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::intersection::~intersection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::intersection::generate() const
{
    labelList nPointIntersections;
    labelList nCellIntersections;
    labelList nFaceIntersections;
    boolList flipMap;

    forAll(zoneGenerators_, i)
    {
        zoneSet zs(zoneGenerators_[i].generate());

        if
        (
            zoneType_ == zoneTypesAll::point
         || (zoneType_ == zoneTypesAll::all && zs.pValid())
        )
        {
            nPointIntersections.setSize(mesh_.nPoints(), 0);
            countIntersections(nPointIntersections, zs.pZone());
        }

        if
        (
            zoneType_ == zoneTypesAll::cell
         || (zoneType_ == zoneTypesAll::all && zs.cValid())
        )
        {
            nCellIntersections.setSize(mesh_.nCells(), 0);
            countIntersections(nCellIntersections, zs.cZone());
        }

        if
        (
            zoneType_ == zoneTypesAll::face
         || (zoneType_ == zoneTypesAll::all && zs.fValid())
        )
        {
            nFaceIntersections.setSize(mesh_.nFaces(), 0);

            if (zs.fZone().oriented())
            {
                flipMap.setSize(mesh_.nFaces(), false);
                countIntersections
                (
                    nFaceIntersections,
                    flipMap,
                    zs.fZone(),
                    zs.fZone().flipMap()
                );
            }
            else
            {
                countIntersections(nFaceIntersections, zs.fZone());
            }
        }
    }

    moveUpdate_ = zoneGenerators_.moveUpdate();

    const labelList faceIndices(indices(nFaceIntersections));

    return zoneSet
    (
        nPointIntersections.size()
      ? new pointZone
        (
            zoneName_,
            indices(nPointIntersections),
            mesh_.pointZones(),
            moveUpdate_,
            true
        )
      : nullptr,
        nCellIntersections.size()
      ? new cellZone
        (
            zoneName_,
            indices(nCellIntersections),
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
      : nullptr,
        nFaceIntersections.size()
      ? flipMap.size()
          ? new faceZone
            (
                zoneName_,
                faceIndices,
                boolList(flipMap, faceIndices),
                mesh_.faceZones(),
                moveUpdate_,
                true
            )
          : new faceZone
            (
                zoneName_,
                faceIndices,
                mesh_.faceZones(),
                moveUpdate_,
                true
            )
      : nullptr
    );
}


// ************************************************************************* //
