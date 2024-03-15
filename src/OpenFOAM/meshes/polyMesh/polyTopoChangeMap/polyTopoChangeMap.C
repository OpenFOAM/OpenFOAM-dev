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
    labelList&& pointMap,
    List<objectMap>&& pointsFromPoints,
    labelList&& faceMap,
    List<objectMap>&& facesFromFaces,
    labelList&& cellMap,
    List<objectMap>&& cellsFromCells,
    labelList&& reversePointMap,
    labelList&& reverseFaceMap,
    labelList&& reverseCellMap,
    labelHashSet&& flipFaceFlux,
    labelListList&& patchPointMap,
    labelList&& oldPatchSizes,
    labelList&& oldPatchStarts,
    labelList&& oldPatchNMeshPoints,
    autoPtr<scalarField>&& oldCellVolumesPtr
)
:
    mesh_(mesh),
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    pointMap_(move(pointMap)),
    pointsFromPointsMap_(move(pointsFromPoints)),
    faceMap_(move(faceMap)),
    facesFromFacesMap_(move(facesFromFaces)),
    cellMap_(move(cellMap)),
    cellsFromCellsMap_(move(cellsFromCells)),
    reversePointMap_(move(reversePointMap)),
    reverseFaceMap_(move(reverseFaceMap)),
    reverseCellMap_(move(reverseCellMap)),
    flipFaceFlux_(move(flipFaceFlux)),
    patchPointMap_(move(patchPointMap)),
    oldPatchSizes_(move(oldPatchSizes)),
    oldPatchStarts_(move(oldPatchStarts)),
    oldPatchNMeshPoints_(move(oldPatchNMeshPoints)),
    oldCellVolumesPtr_(move(oldCellVolumesPtr))
{}


// ************************************************************************* //
