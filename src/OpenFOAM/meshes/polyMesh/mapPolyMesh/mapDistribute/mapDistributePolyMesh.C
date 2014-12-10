/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "mapDistributePolyMesh.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::calcPatchSizes()
{
    oldPatchSizes_.setSize(oldPatchStarts_.size());

    // Calculate old patch sizes
    for (label patchI = 0; patchI < oldPatchStarts_.size() - 1; patchI++)
    {
        oldPatchSizes_[patchI] =
            oldPatchStarts_[patchI + 1] - oldPatchStarts_[patchI];
    }

    // Set the last one by hand
    const label lastPatchID = oldPatchStarts_.size() - 1;

    oldPatchSizes_[lastPatchID] = nOldFaces_ - oldPatchStarts_[lastPatchID];

    if (min(oldPatchSizes_) < 0)
    {
        FatalErrorIn("mapDistributePolyMesh::calcPatchSizes()")
            << "Calculated negative old patch size:" << oldPatchSizes_ << nl
            << "Error in mapping data" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    const polyMesh& mesh,

    // mesh before changes
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const Xfer<labelList>& oldPatchStarts,
    const Xfer<labelList>& oldPatchNMeshPoints,

    // how to subset pieces of mesh to send across
    const Xfer<labelListList>& subPointMap,
    const Xfer<labelListList>& subFaceMap,
    const Xfer<labelListList>& subCellMap,
    const Xfer<labelListList>& subPatchMap,

    // how to reconstruct received mesh
    const Xfer<labelListList>& constructPointMap,
    const Xfer<labelListList>& constructFaceMap,
    const Xfer<labelListList>& constructCellMap,
    const Xfer<labelListList>& constructPatchMap
)
:
    mesh_(mesh),
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    oldPatchSizes_(oldPatchStarts().size()),
    oldPatchStarts_(oldPatchStarts),
    oldPatchNMeshPoints_(oldPatchNMeshPoints),
    pointMap_(mesh.nPoints(), subPointMap, constructPointMap),
    faceMap_(mesh.nFaces(), subFaceMap, constructFaceMap),
    cellMap_(mesh.nCells(), subCellMap, constructCellMap),
    patchMap_(mesh.boundaryMesh().size(), subPatchMap, constructPatchMap)
{
    calcPatchSizes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::distributePointIndices(labelList& lst) const
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


void Foam::mapDistributePolyMesh::distributeFaceIndices(labelList& lst) const
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


void Foam::mapDistributePolyMesh::distributeCellIndices(labelList& lst) const
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


void Foam::mapDistributePolyMesh::distributePatchIndices(labelList& lst) const
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


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
