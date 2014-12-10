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

#include "starMesh.H"
#include "emptyPolyPatch.H"
#include "demandDrivenData.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Merge tolerances
// Moni, tolerances:
const scalar starMesh::smallMergeTol_ = 1e-3;
const scalar starMesh::cpMergePointTol_ = 1e-4;

// const scalar starMesh::smallMergeTol_ = 4e-4;
// const scalar starMesh::cpMergePointTol_ = 1e-3;

// Cell shape models
const cellModel* starMesh::unknownPtr_ = cellModeller::lookup("unknown");
const cellModel* starMesh::tetPtr_ = cellModeller::lookup("tet");
const cellModel* starMesh::pyrPtr_ = cellModeller::lookup("pyr");
const cellModel* starMesh::tetWedgePtr_ = cellModeller::lookup("tetWedge");
const cellModel* starMesh::prismPtr_ = cellModeller::lookup("prism");
const cellModel* starMesh::wedgePtr_ = cellModeller::lookup("wedge");
const cellModel* starMesh::hexPtr_ = cellModeller::lookup("hex");

const cellModel* starMesh::sammTrim1Ptr_ = cellModeller::lookup("sammTrim1");
const cellModel* starMesh::sammTrim2Ptr_ = cellModeller::lookup("sammTrim2");
const cellModel* starMesh::sammTrim3Ptr_ = cellModeller::lookup("sammTrim3");
const cellModel* starMesh::sammTrim4Ptr_ = cellModeller::lookup("sammTrim4");
const cellModel* starMesh::sammTrim5Ptr_ = cellModeller::lookup("sammTrim5");
const cellModel* starMesh::sammTrim8Ptr_ =
    cellModeller::lookup("hexagonalPrism");

// Regular cell point addressing
// SAMM point addressing
const label starMesh::regularAddressingTable[6][8] =
{
    { 0,  1,  2,  4, -1, -1, -1, -1},    // tet
    { 0,  1,  2,  3,  4, -1, -1, -1},    // pyramid
    { 0,  1,  2,  4,  6, -1, -1, -1},    // tet wedge
    { 0,  1,  2,  4,  5,  6, -1, -1},    // prism
    { 7,  6,  5,  3,  2,  1,  0, -1},    // wedge
    { 0,  1,  2,  3,  4,  5,  6,  7}     // hex
};


// SAMM point addressing
const label starMesh::sammAddressingTable[9][12] =
{
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},    // samm0 - empty
    { 3,  2,  6,  7, 11,  9,  1,  5,  4, 12, -1, -1},    // samm1+
    {13,  5,  6,  2, 10, 12,  4,  7,  3, 11, -1, -1},    // samm2+
    { 2,  3,  0,  1, 10, 11, 12,  4,  8,  9, -1, -1},    // samm3+
    { 0,  1,  3,  4, 13,  8,  9, 10, 11, 12, -1, -1},    // samm4+
    {12,  7,  6,  5,  8, 11, 10,  9, -1, -1, -1, -1},    // samm5+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},    // samm6 - empty
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},    // samm7 - empty
    {11,  3, 15, 12,  4,  8, 10,  2, 14, 13,  5,  9}     // samm8+
};


// lookup table giving OpenFOAM face number when looked up with shape index
// (first index) and STAR face number
// - first column is always -1
// - last column is -1 for all but hexagonal prism
// WARNING: Possible bug for sammTrim2
// The lookup table for SAMM shapes is based on the rotation of the
// shape. This would imply that the table below needs to be split between
// the regular shapes (3-9), which are OK, and the SAMM shapes, for which
// the face lookup needs to be done based on the rotation. Thus, for a samm
// cell, firts find out the face index in the normal rotation using the cell
// face permutation table and then use the index from the shape face lookup.
// Additionally, have in mind that this silliness does not allow matches
// on face 7 and 8 of the samm cell.

const label starMesh::sammFacePermutationTable[24][8] =
{
  {-1, 1, 2, 3, 4, 5, 6, 7},    // permutation   0
  {-1, 3, 4, 5, 6, 1, 2, 7},    // permutation   1
  {-1, 5, 6, 1, 2, 3, 4, 7},    // permutation   2
  {-1, 1, 2, 5, 6, 4, 3, 7},    // permutation   3
  {-1, 3, 4, 1, 2, 6, 5, 7},    // permutation   4
  {-1, 5, 6, 3, 4, 2, 1, 7},    // permutation   5
  {-1, 1, 2, 4, 3, 6, 5, 7},    // permutation   6
  {-1, 3, 4, 6, 5, 2, 1, 7},    // permutation   7
  {-1, 5, 6, 2, 1, 4, 3, 7},    // permutation   8
  {-1, 1, 2, 6, 5, 3, 4, 7},    // permutation   9
  {-1, 3, 4, 2, 1, 5, 6, 7},    // permutation  10
  {-1, 5, 6, 4, 3, 1, 2, 7},    // permutation  11
  {-1, 2, 1, 5, 6, 3, 4, 7},    // permutation  12
  {-1, 4, 3, 1, 2, 5, 6, 7},    // permutation  13
  {-1, 6, 5, 3, 4, 1, 2, 7},    // permutation  14
  {-1, 2, 1, 3, 4, 6, 5, 7},    // permutation  15
  {-1, 4, 3, 5, 6, 2, 1, 7},    // permutation  16
  {-1, 6, 5, 1, 2, 4, 3, 7},    // permutation  17
  {-1, 2, 1, 6, 5, 4, 3, 7},    // permutation  18
  {-1, 4, 3, 2, 1, 6, 5, 7},    // permutation  19
  {-1, 6, 5, 4, 3, 2, 1, 7},    // permutation  20
  {-1, 2, 1, 4, 3, 5, 6, 7},    // permutation  21
  {-1, 4, 3, 6, 5, 1, 2, 7},    // permutation  22
  {-1, 6, 5, 2, 1, 3, 4, 7}     // permutation  23
};

const label starMesh::shapeFaceLookup[19][9] =
{
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  0 - empty+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  1 - empty+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  2 - empty+
    {-1,  4,  5,  2,  3,  0,  1, -1, -1},    // shape  3 - hex+
    {-1,  1,  0,  5,  4,  2,  3, -1, -1},    // shape  4 - wedge+
    {-1,  0,  1,  4, -1,  2,  3, -1, -1},    // shape  5 - prism+
    {-1,  0, -1,  4,  2,  1,  3, -1, -1},    // shape  6 - pyr+
    {-1,  3, -1,  2, -1,  1,  0, -1, -1},    // shape  7 - tet+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  8 - splitHex (empty)
    {-1,  0, -1,  1, -1,  2,  3, -1, -1},    // shape  9 - tetWedge+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape 10 - empty+
    {-1,  1,  0,  3,  2,  5,  4,  6, -1},    // shape 11 - sammTrim1
    {-1,  5,  4,  1,  0,  3,  2,  6, -1},    // shape 12 - sammTrim2
    {-1,  2,  3,  0,  1,  4,  5,  6, -1},    // shape 13 - sammTrim3
    {-1,  2,  3,  0,  1,  4,  5,  6, -1},    // shape 14 - sammTrim4
    {-1,  5,  2,  4,  3,  1,  0, -1, -1},    // shape 15 - sammTrim5
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape 16 - empty
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape 17 - empty
    {-1,  1,  0,  6,  7,  2,  3,  4,  5}     // shape 18 - sammTrim8
};


// The star to foam face order mapping tables are potentially incomplete
// Currently available data is listed below.
// 1) hex and degenerate hex: OK
// samm trim 1:
// star number: 1 2 3 4 5 6 7 8  In ROTATION 0
// foam number: 5 4 1 0 3 2 6
// confirmed:   1 0 3 2 5 4 6

// samm trim 2:
// star number: 1 2 3 4 5 6 7 8  In ROTATION 0
// foam number: 5 4 1 0 3 2 6
// confirmed:     4   0 3 2

// samm trim 3:
// star number: 1 2 3 4 5 6 7 8  In ROTATION 0
// foam number: 2 3 0 1 4 5 6
// confirmed:

// samm trim 4:
// star number: 1 2 3 4 5 6 7 8  In ROTATION 0
// foam number: 2 3 0 1 4 5 6
// confirmed:

// samm trim 5:
// star number: 1 2 3 4 5 6 7 8  In ROTATION 0
// foam number: 2 4 3 1 0 5
// confirmed:   5 2 4 3 1 0

// samm trim 8:
// star number: 1 2 3 4 5 6 7 8  In ROTATION 0
// foam number: 2 5 4 7 1 0 3 6
// confirmed:   1 0 6


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Make polyhedral mesh data (packing)
void starMesh::createPolyMeshData()
{
    Info<< "Creating a polyMesh" << endl;

    createPolyCells();

    Info<< "\nNumber of internal faces: "
        << nInternalFaces_ << endl;

    createPolyBoundary();
}


// Clear extra storage before creation of the mesh to remove
// a memory peak
void starMesh::clearExtraStorage()
{
    Info<< "Clearing extra storage" << endl;

    starPointLabelLookup_.setSize(0);
    starPointID_.setSize(0);
    starCellID_.setSize(0);
    starCellLabelLookup_.setSize(0);
    starCellPermutation_.setSize(0);
    cellFaces_.setSize(0);
    boundaryCellIDs_.setSize(0);
    boundaryCellFaceIDs_.setSize(0);
    couples_.clear();

    deleteDemandDrivenData(pointCellsPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
starMesh::starMesh
(
    const fileName& prefix,
    const Time& rt,
    const scalar scaleFactor
)
:
    casePrefix_(prefix),
    runTime_(rt),
    points_(0),
    cellShapes_(0),
    boundary_(0),
    patchTypes_(0),
    defaultFacesName_("defaultFaces"),
    defaultFacesType_(emptyPolyPatch::typeName),
    patchNames_(0),
    patchPhysicalTypes_(0),
    starPointLabelLookup_(0),
    starPointID_(0),
    starCellID_(0),
    starCellLabelLookup_(0),
    starCellPermutation_(0),
    cellFaces_(0),
    boundaryCellIDs_(0),
    boundaryCellFaceIDs_(0),
    meshFaces_(0),
    cellPolys_(0),
    nInternalFaces_(0),
    polyBoundaryPatchStartIndices_(0),
    pointCellsPtr_(NULL),
    couples_(0),
    isShapeMesh_(true)
{
    readPoints(scaleFactor);

    readCells();

    readBoundary();

    fixCollapsedEdges();

    readCouples();

    if (couples_.size())
    {
        createCoupleMatches();
    }

    markBoundaryFaces();

    mergeCoupleFacePoints();

    purgeCellShapes();

    collectBoundaryFaces();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

starMesh::~starMesh()
{
    deleteDemandDrivenData(pointCellsPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
