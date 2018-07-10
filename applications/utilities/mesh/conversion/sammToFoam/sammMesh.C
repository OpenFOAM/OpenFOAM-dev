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

#include "sammMesh.H"
#include "emptyPolyPatch.H"
#include "demandDrivenData.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Cell shape models
const Foam::cellModel* Foam::sammMesh::unknownPtr_ =
    Foam::cellModeller::lookup("unknown");
const Foam::cellModel* Foam::sammMesh::hexPtr_ =
    Foam::cellModeller::lookup("hex");
const Foam::cellModel* Foam::sammMesh::wedgePtr_ =
    Foam::cellModeller::lookup("wedge");
const Foam::cellModel* Foam::sammMesh::prismPtr_ =
    Foam::cellModeller::lookup("prism");
const Foam::cellModel* Foam::sammMesh::pyrPtr_ =
    Foam::cellModeller::lookup("pyr");
const Foam::cellModel* Foam::sammMesh::tetPtr_ =
    Foam::cellModeller::lookup("tet");
const Foam::cellModel* Foam::sammMesh::tetWedgePtr_ =
    Foam::cellModeller::lookup("tetWedge");

const Foam::cellModel* Foam::sammMesh::sammTrim1Ptr_ =
    Foam::cellModeller::lookup("sammTrim1");
const Foam::cellModel* Foam::sammMesh::sammTrim2Ptr_ =
    Foam::cellModeller::lookup("sammTrim2");
const Foam::cellModel* Foam::sammMesh::sammTrim3Ptr_ =
    Foam::cellModeller::lookup("sammTrim3");
const Foam::cellModel* Foam::sammMesh::sammTrim4Ptr_ =
    Foam::cellModeller::lookup("sammTrim4");
const Foam::cellModel* Foam::sammMesh::sammTrim5Ptr_ =
    Foam::cellModeller::lookup("sammTrim5");
const Foam::cellModel* Foam::sammMesh::sammTrim8Ptr_ =
    Foam::cellModeller::lookup("hexagonalPrism");

// lookup table giving OpenFOAM face number when looked up with shape index
// (first index) and STAR face number
// - first column is always -1
// - last column is -1 for all but hexagonal prism
// WARNING: Possible bug for sammTrim2
// There is a possibility that the lookup table for SAMM shapes is based on
// the rotation of the shape. This would imply that the table below would need
// to be split between the regular shapes (3-9), which are OK, and the SAMM
// shapes, for which the face lookup needs to be done based on the rotation.
// However, at the moment I haven't got enough info to complete the toble and
// there are no cases that break it. Please reconsider in the light of mode
// information.
const Foam::label Foam::sammMesh::shapeFaceLookup[19][9] =
{
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  0 - empty+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  1 - empty+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  2 - empty+
    {-1,  4,  5,  2,  3,  0,  1, -1, -1},    // shape  3 - hex+
    {-1,  4,  5,  2,  3,  0,  1, -1, -1},    // shape  4 - wedge+
    {-1,  0,  1,  4, -1,  2,  3, -1, -1},    // shape  5 - prism+
    {-1,  0, -1,  4,  2,  1,  3, -1, -1},    // shape  6 - pyr+
    {-1,  3, -1,  2, -1,  1,  0, -1, -1},    // shape  7 - tet+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape  8 - splitHex (empty)
    {-1,  0, -1,  1, -1,  2,  3, -1, -1},    // shape  9 - tetWedge+
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape 10 - empty+
    {-1,  5,  4,  0,  1,  2,  3,  6, -1},    // shape 11 - sammTrim1+
//    {-1,  1,  0,  2,  3,  4,  5,  6, -1},    // shape 12 - sammTrim2 ?
    {-1, 1,  0,  2,  4,  3,  5,  6, -1},    // shape 12 - sammTrim2  f(4)=4
    {-1,  5,  4,  0,  1,  2,  3,  6, -1},    // shape 13 - sammTrim3+
    {-1,  5,  4,  1,  0,  3,  2,  6, -1},    // shape 14 - sammTrim4
    {-1,  4,  3,  2,  5,  1,  0, -1, -1},    // shape 15 - sammTrim5
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape 16 - empty
    {-1, -1, -1, -1, -1, -1, -1, -1, -1},    // shape 17 - empty
    {-1,  0,  1,  2,  5,  3,  6,  4,  7}     // shape 18 - sammTrim8
};

// SAMM cell lookup data

// List of pointers used instead of pointer list o avoid
// de-allocation problems
Foam::List<const Foam::cellModel*> Foam::sammMesh::sammShapeLookup
(
    256,
    reinterpret_cast<cellModel*>(0)
);

Foam::List<const Foam::label*> Foam::sammMesh::sammAddressingTable
(
    256,
    reinterpret_cast<label*>(0)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sammMesh::createPolyMeshData()
{
    Info<< "Creating a polyMesh" << endl;

    createPolyCells();

    Info<< "\nNumber of internal faces: "
        << nInternalFaces_ << endl;

    createPolyBoundary();

    label nProblemCells = 0;

    // check that there is no zeros in the cellPolys_
    forAll(cellPolys_, celli)
    {
        const labelList& curFaceLabels = cellPolys_[celli];

        forAll(curFaceLabels, facei)
        {
            if (curFaceLabels[facei] == -1)
            {
                Info<< "cell " << celli
                    << " has got an unmatched face. "
                    << "Index: " << cellShapes_[celli].model().index() << endl
//                     << "cell shape: " << cellShapes_[celli] << endl
//                     << "shape faces: " << cellShapes_[celli].faces() << endl
                    << "cellPolys: " << cellPolys_[celli] << endl
//                     << "cell faces: " << cellFaces_[celli]
                    << endl;

                nProblemCells++;

                break;
            }
        }
    }

    if (nProblemCells > 0)
    {
        Info<< "Number of problem cells: " << nProblemCells << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sammMesh::sammMesh
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
    starCellLabelLookup_(0),
    cellFaces_(0),
    meshFaces_(0),
    cellPolys_(0),
    nInternalFaces_(0),
    polyBoundaryPatchStartIndices_(0),
    pointCellsPtr_(nullptr),
    isShapeMesh_(true)
{
    // Fill in the lookup tables
    fillSammCellShapeTable();
    fillSammAddressingTable();

    readPoints(scaleFactor);

    readCells();

    readBoundary();

    fixCollapsedEdges();

    readCouples();

    // create boundary faces
    createBoundaryFaces();

    // after all this is done do couples
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sammMesh::~sammMesh()
{
    deleteDemandDrivenData(pointCellsPtr_);
}


// ************************************************************************* //
