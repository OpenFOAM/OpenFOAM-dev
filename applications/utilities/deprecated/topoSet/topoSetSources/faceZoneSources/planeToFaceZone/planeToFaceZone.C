/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "planeToFaceZone.H"
#include "polyMesh.H"
#include "patchRegionSplit.H"
#include "faceZoneSet.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::NamedEnum<Foam::planeToFaceZone::include, 2>
Foam::planeToFaceZone::includeNames_
{
    "all",
    "closest"
};

namespace Foam
{
    defineTypeNameAndDebug(planeToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, planeToFaceZone, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::planeToFaceZone::combine(faceZoneSet& fzSet, const bool add) const
{
    // Mark all cells with centres above the plane
    boolList cellIsAbovePlane(mesh_.nCells());
    forAll(mesh_.cells(), celli)
    {
        cellIsAbovePlane[celli] =
            ((mesh_.cellCentres()[celli] - point_) & normal_) > 0;
    }

    // Mark all coupled neighbour cells with centres above the plane
    boolList bFaceNbrCellIsAbovePlane(mesh_.nFaces() - mesh_.nInternalFaces());
    {
        vectorField bFaceNbrCellCentres;
        syncTools::swapBoundaryCellPositions
        (
            mesh_,
            mesh_.cellCentres(),
            bFaceNbrCellCentres
        );
        forAll(bFaceNbrCellIsAbovePlane, bFacei)
        {
            bFaceNbrCellIsAbovePlane[bFacei] =
                ((bFaceNbrCellCentres[bFacei] - point_) & normal_) > 0;
        }
    }

    // Mark all faces that sit between cells above and below the plane
    boolList faceIsOnPlane(mesh_.nFaces(), false);
    forAll(mesh_.faceNeighbour(), facei)
    {
        faceIsOnPlane[facei] =
            cellIsAbovePlane[mesh_.faceOwner()[facei]]
         != cellIsAbovePlane[mesh_.faceNeighbour()[facei]];
    }
    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& patch = mesh_.boundary()[patchi];
        if (!patch.coupled()) continue;
        forAll(patch, patchFacei)
        {
            const label facei = patch.start() + patchFacei;
            faceIsOnPlane[facei] =
                cellIsAbovePlane[mesh_.faceOwner()[facei]]
             != bFaceNbrCellIsAbovePlane[facei - mesh_.nInternalFaces()];
        }
    }

    // Ensure consistency across couplings
    syncTools::syncFaceList(mesh_, faceIsOnPlane, orEqOp());

    // Convert marked faces to a list of indices
    labelList newSetFaces(findIndices(faceIsOnPlane, true));

    // If constructing a single contiguous set, remove all faces except those
    // connected to the contiguous region closest to the specified point
    if (include_ == include::closest)
    {
        // Identify contiguous regions
        const patchRegionSplit prs(mesh_, newSetFaces);

        // Find the smallest distance from each region to the point
        scalarField regionMinDistSqr(prs.nRegions(), vGreat);
        forAll(newSetFaces, fi)
        {
            const label facei = newSetFaces[fi];
            const label regioni = prs[fi];

            const scalar distSqr =
                magSqr
                (
                    mesh_.faces()[facei].nearestPoint
                    (
                        point_,
                        mesh_.points()
                    ).rawPoint()
                  - point_
                );

            regionMinDistSqr[regioni] = min(regionMinDistSqr[regioni], distSqr);
        }
        Pstream::listCombineGather(regionMinDistSqr, minEqOp());
        Pstream::listCombineScatter(regionMinDistSqr);

        // Choose the closest region
        const label selectedRegioni = findMin(regionMinDistSqr);

        // Remove faces in other regions from the list by shuffling up
        label fi0 = 0;
        forAll(newSetFaces, fi)
        {
            newSetFaces[fi0] = newSetFaces[fi];
            if (prs[fi] == selectedRegioni) fi0 ++;
        }
        newSetFaces.resize(fi0);
    }

    // Modify the face zone set
    DynamicList<label> newAddressing;
    DynamicList<bool> newFlipMap;
    if (add)
    {
        // Start from copy
        newAddressing = DynamicList<label>(fzSet.addressing());
        newFlipMap = DynamicList<bool>(fzSet.flipMap());

        // Add anything from the new set that is not already in the zone set
        forAll(newSetFaces, newSetFacei)
        {
            const label facei = newSetFaces[newSetFacei];

            if (!fzSet.found(facei))
            {
                newAddressing.append(facei);
                newFlipMap.append(cellIsAbovePlane[mesh_.faceOwner()[facei]]);
            }
        }
    }
    else
    {
        // Start from empty
        newAddressing = DynamicList<label>(fzSet.addressing().size());
        newFlipMap = DynamicList<bool>(fzSet.flipMap().size());

        // Add everything from the zone set that is not also in the new set
        labelHashSet newSet(newSetFaces);
        forAll(fzSet.addressing(), i)
        {
            const label facei = fzSet.addressing()[i];

            if (!newSet.found(facei))
            {
                newAddressing.append(facei);
                newFlipMap.append(cellIsAbovePlane[mesh_.faceOwner()[facei]]);
            }
        }
    }
    fzSet.addressing().transfer(newAddressing);
    fzSet.flipMap().transfer(newFlipMap);
    fzSet.updateSet();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::planeToFaceZone::planeToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    point_(dict.lookup<vector>("point", dimensions::length)),
    normal_(dict.lookup<vector>("normal", dimless)),
    include_
    (
        includeNames_
        [
            dict.lookupOrDefault<word>
            (
                "include",
                includeNames_[include::all]
            )
        ]
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::planeToFaceZone::~planeToFaceZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::planeToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet." << endl;
    }
    else
    {
        faceZoneSet& fzSet = refCast<faceZoneSet>(set);

        if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
        {
            Info<< "    Adding faces which form a plane at " << point_
                << " with normal " << normal_ << endl;

            combine(fzSet, true);
        }
        else if (action == topoSetSource::DELETE)
        {
            Info<< "    Removing faces which form a plane at " << point_
                << " with normal " << normal_ << endl;

            combine(fzSet, false);
        }
    }
}


// ************************************************************************* //
