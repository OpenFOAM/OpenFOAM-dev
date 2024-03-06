/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "polyTopoChangeMap.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyTopoChangeMap::polyTopoChangeMap(const polyMesh& mesh)
:
    mesh_(mesh),
    nOldPoints_(-1),
    nOldFaces_(-1),
    nOldCells_(-1)
{}


Foam::polyTopoChangeMap::polyTopoChangeMap
(
    const polyMesh& mesh,
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const labelList& pointMap,
    const List<objectMap>& pointsFromPoints,
    const labelList& faceMap,
    const List<objectMap>& facesFromPoints,
    const List<objectMap>& facesFromEdges,
    const List<objectMap>& facesFromFaces,
    const labelList& cellMap,
    const List<objectMap>& cellsFromCells,
    const labelList& reversePointMap,
    const labelList& reverseFaceMap,
    const labelList& reverseCellMap,
    const labelHashSet& flipFaceFlux,
    const labelListList& patchPointMap,
    const pointField& preMotionPoints,
    const labelList& oldPatchStarts,
    const labelList& oldPatchNMeshPoints,
    const autoPtr<scalarField>& oldCellVolumesPtr
)
:
    mesh_(mesh),
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    pointMap_(pointMap),
    pointsFromPointsMap_(pointsFromPoints),
    faceMap_(faceMap),
    facesFromPointsMap_(facesFromPoints),
    facesFromEdgesMap_(facesFromEdges),
    facesFromFacesMap_(facesFromFaces),
    cellMap_(cellMap),
    cellsFromCellsMap_(cellsFromCells),
    reversePointMap_(reversePointMap),
    reverseFaceMap_(reverseFaceMap),
    reverseCellMap_(reverseCellMap),
    flipFaceFlux_(flipFaceFlux),
    patchPointMap_(patchPointMap),
    preMotionPoints_(preMotionPoints),
    oldPatchSizes_(oldPatchStarts.size()),
    oldPatchStarts_(oldPatchStarts),
    oldPatchNMeshPoints_(oldPatchNMeshPoints),
    oldCellVolumesPtr_(oldCellVolumesPtr)
{
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

        if (polyMesh::debug)
        {
            if (min(oldPatchSizes_) < 0)
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }
        }
    }
}


Foam::polyTopoChangeMap::polyTopoChangeMap
(
    const polyMesh& mesh,
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    labelList& pointMap,
    List<objectMap>& pointsFromPoints,
    labelList& faceMap,
    List<objectMap>& facesFromPoints,
    List<objectMap>& facesFromEdges,
    List<objectMap>& facesFromFaces,
    labelList& cellMap,
    List<objectMap>& cellsFromCells,
    labelList& reversePointMap,
    labelList& reverseFaceMap,
    labelList& reverseCellMap,
    labelHashSet& flipFaceFlux,
    labelListList& patchPointMap,
    pointField& preMotionPoints,
    labelList& oldPatchStarts,
    labelList& oldPatchNMeshPoints,
    autoPtr<scalarField>& oldCellVolumesPtr,
    const bool reuse
)
:
    mesh_(mesh),
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    pointMap_(pointMap, reuse),
    pointsFromPointsMap_(pointsFromPoints, reuse),
    faceMap_(faceMap, reuse),
    facesFromPointsMap_(facesFromPoints, reuse),
    facesFromEdgesMap_(facesFromEdges, reuse),
    facesFromFacesMap_(facesFromFaces, reuse),
    cellMap_(cellMap, reuse),
    cellsFromCellsMap_(cellsFromCells, reuse),
    reversePointMap_(reversePointMap, reuse),
    reverseFaceMap_(reverseFaceMap, reuse),
    reverseCellMap_(reverseCellMap, reuse),
    flipFaceFlux_(flipFaceFlux),
    patchPointMap_(patchPointMap, reuse),
    preMotionPoints_(preMotionPoints, reuse),
    oldPatchSizes_(oldPatchStarts.size()),
    oldPatchStarts_(oldPatchStarts, reuse),
    oldPatchNMeshPoints_(oldPatchNMeshPoints, reuse),
    oldCellVolumesPtr_(oldCellVolumesPtr, reuse)
{
    if (oldPatchStarts_.size() > 0)
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

        if (polyMesh::debug)
        {
            if (min(oldPatchSizes_) < 0)
            {
                FatalErrorInFunction
                    << "Calculated negative old patch size."
                    << "  Error in mapping data"
                    << abort(FatalError);
            }
        }
    }
}


// ************************************************************************* //
