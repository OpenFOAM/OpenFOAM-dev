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

#include "polyDistributionMap.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyDistributionMap::calcPatchSizes()
{
    oldPatchSizes_.setSize(oldPatchStarts_.size());

    if (oldPatchStarts_.size())
    {
        // Calculate old patch sizes
        for (label patchi = 0; patchi < oldPatchStarts_.size() - 1; patchi++)
        {
            oldPatchSizes_[patchi] =
                oldPatchStarts_[patchi + 1] - oldPatchStarts_[patchi];
        }

        // Set the last one by hand
        const label lastPatchID = oldPatchStarts_.size() - 1;

        oldPatchSizes_[lastPatchID] = nOldFaces_ - oldPatchStarts_[lastPatchID];

        if (min(oldPatchSizes_) < 0)
        {
            FatalErrorInFunction
                << "Calculated negative old patch size:" << oldPatchSizes_ << nl
                << "Error in mapping data" << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyDistributionMap::polyDistributionMap
(
    const polyMesh& mesh,

    // mesh before changes
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    labelList&& oldPatchStarts,
    labelList&& oldPatchNMeshPoints,

    // how to subset pieces of mesh to send across
    labelListList&& subPointMap,
    labelListList&& subFaceMap,
    labelListList&& subCellMap,
    labelListList&& subPatchMap,

    // how to reconstruct received mesh
    labelListList&& constructPointMap,
    labelListList&& constructFaceMap,
    labelListList&& constructCellMap,
    labelListList&& constructPatchMap,

    const bool subFaceHasFlip,
    const bool constructFaceHasFlip
)
:
    mesh_(mesh),
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    oldPatchSizes_(oldPatchStarts.size()),
    oldPatchStarts_(move(oldPatchStarts)),
    oldPatchNMeshPoints_(move(oldPatchNMeshPoints)),
    pointMap_(mesh.nPoints(), move(subPointMap), move(constructPointMap)),
    faceMap_
    (
        mesh.nFaces(),
        move(subFaceMap),
        move(constructFaceMap),
        move(subFaceHasFlip),
        constructFaceHasFlip
    ),
    cellMap_(mesh.nCells(), move(subCellMap), move(constructCellMap)),
    patchMap_
    (
        mesh.boundaryMesh().size(),
        move(subPatchMap),
        move(constructPatchMap)
    )
{
    calcPatchSizes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyDistributionMap::distributePointIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        createWithValues<boolList>
        (
            nOldPoints(),
            false,
            lst,
            true
        )
    );

    // Distribute
    distributePointData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


void Foam::polyDistributionMap::distributeFaceIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        createWithValues<boolList>
        (
            nOldFaces(),
            false,
            lst,
            true
        )
    );

    // Distribute
    distributeFaceData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


void Foam::polyDistributionMap::distributeCellIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        createWithValues<boolList>
        (
            nOldCells(),
            false,
            lst,
            true
        )
    );

    // Distribute
    distributeCellData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


void Foam::polyDistributionMap::distributePatchIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        createWithValues<boolList>
        (
            oldPatchStarts().size(),    // nOldPatches
            false,
            lst,
            true
        )
    );

    // Distribute
    distributePatchData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


// ************************************************************************* //
