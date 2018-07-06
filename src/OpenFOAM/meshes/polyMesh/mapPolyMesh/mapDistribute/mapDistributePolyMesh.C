/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "mapDistributePolyMesh.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::calcPatchSizes()
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

Foam::mapDistributePolyMesh::mapDistributePolyMesh()
:
    nOldPoints_(0),
    nOldFaces_(0),
    nOldCells_(0),
    oldPatchSizes_(0),
    oldPatchStarts_(0),
    oldPatchNMeshPoints_(0),
    pointMap_(),
    faceMap_(),
    cellMap_(),
    patchMap_()
{}


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
    const Xfer<labelListList>& constructPatchMap,

    const bool subFaceHasFlip,
    const bool constructFaceHasFlip
)
:
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    oldPatchSizes_(oldPatchStarts().size()),
    oldPatchStarts_(oldPatchStarts),
    oldPatchNMeshPoints_(oldPatchNMeshPoints),
    pointMap_(mesh.nPoints(), subPointMap, constructPointMap),
    faceMap_
    (
        mesh.nFaces(),
        subFaceMap,
        constructFaceMap,
        subFaceHasFlip,
        constructFaceHasFlip
    ),
    cellMap_(mesh.nCells(), subCellMap, constructCellMap),
    patchMap_(mesh.boundaryMesh().size(), subPatchMap, constructPatchMap)
{
    calcPatchSizes();
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    // mesh before changes
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const Xfer<labelList>& oldPatchStarts,
    const Xfer<labelList>& oldPatchNMeshPoints,

    // how to transfer pieces of mesh
    const Xfer<mapDistribute>& pointMap,
    const Xfer<mapDistribute>& faceMap,
    const Xfer<mapDistribute>& cellMap,
    const Xfer<mapDistribute>& patchMap
)
:
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    oldPatchSizes_(oldPatchStarts().size()),
    oldPatchStarts_(oldPatchStarts),
    oldPatchNMeshPoints_(oldPatchNMeshPoints),
    pointMap_(pointMap),
    faceMap_(faceMap),
    cellMap_(cellMap),
    patchMap_(patchMap)
{
    calcPatchSizes();
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    const Xfer<mapDistributePolyMesh>& map
)
:
    nOldPoints_(map().nOldPoints_),
    nOldFaces_(map().nOldFaces_),
    nOldCells_(map().nOldCells_),
    oldPatchSizes_(map().oldPatchSizes_.xfer()),
    oldPatchStarts_(map().oldPatchStarts_.xfer()),
    oldPatchNMeshPoints_(map().oldPatchNMeshPoints_.xfer()),
    pointMap_(map().pointMap_.xfer()),
    faceMap_(map().faceMap_.xfer()),
    cellMap_(map().cellMap_.xfer()),
    patchMap_(map().patchMap_.xfer())
{}


Foam::mapDistributePolyMesh::mapDistributePolyMesh(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::transfer(mapDistributePolyMesh& rhs)
{
    nOldPoints_ = rhs.nOldPoints_;
    nOldFaces_ = rhs.nOldFaces_;
    nOldCells_ = rhs.nOldCells_;
    oldPatchSizes_.transfer(rhs.oldPatchSizes_);
    oldPatchStarts_.transfer(rhs.oldPatchStarts_);
    oldPatchNMeshPoints_.transfer(rhs.oldPatchNMeshPoints_);
    pointMap_.transfer(rhs.pointMap_);
    faceMap_.transfer(rhs.faceMap_);
    cellMap_.transfer(rhs.cellMap_);
    patchMap_.transfer(rhs.patchMap_);
}


Foam::Xfer<Foam::mapDistributePolyMesh> Foam::mapDistributePolyMesh::xfer()
{
    return xferMove(*this);
}


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


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::operator=(const mapDistributePolyMesh& rhs)
{
    nOldPoints_ = rhs.nOldPoints_;
    nOldFaces_ = rhs.nOldFaces_;
    nOldCells_ = rhs.nOldCells_;
    oldPatchSizes_ = rhs.oldPatchSizes_;
    oldPatchStarts_ = rhs.oldPatchStarts_;
    oldPatchNMeshPoints_ = rhs.oldPatchNMeshPoints_;
    pointMap_ = rhs.pointMap_;
    faceMap_ = rhs.faceMap_;
    cellMap_ = rhs.cellMap_;
    patchMap_ = rhs.patchMap_;
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, mapDistributePolyMesh& map)
{
    is.fatalCheck("operator>>(Istream&, mapDistributePolyMesh&)");

    is  >> map.nOldPoints_
        >> map.nOldFaces_
        >> map.nOldCells_
        >> map.oldPatchSizes_
        >> map.oldPatchStarts_
        >> map.oldPatchNMeshPoints_
        >> map.pointMap_
        >> map.faceMap_
        >> map.cellMap_
        >> map.patchMap_;

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const mapDistributePolyMesh& map)
{
    os  << map.nOldPoints_
        << token::SPACE << map.nOldFaces_
        << token::SPACE << map.nOldCells_ << token::NL
        << map.oldPatchSizes_ << token::NL
        << map.oldPatchStarts_ << token::NL
        << map.oldPatchNMeshPoints_ << token::NL
        << map.pointMap_ << token::NL
        << map.faceMap_ << token::NL
        << map.cellMap_ << token::NL
        << map.patchMap_;

    return os;
}


// ************************************************************************* //
