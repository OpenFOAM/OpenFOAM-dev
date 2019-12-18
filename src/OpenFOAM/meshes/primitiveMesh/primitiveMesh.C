/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(primitiveMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primitiveMesh::primitiveMesh()
:
    nInternalPoints_(0),    // note: points are considered ordered on empty mesh
    nPoints_(0),
    nInternal0Edges_(-1),
    nInternal1Edges_(-1),
    nInternalEdges_(-1),
    nEdges_(-1),
    nInternalFaces_(0),
    nFaces_(0),
    nCells_(0),

    cellShapesPtr_(nullptr),
    edgesPtr_(nullptr),
    ccPtr_(nullptr),
    ecPtr_(nullptr),
    pcPtr_(nullptr),

    cfPtr_(nullptr),
    efPtr_(nullptr),
    pfPtr_(nullptr),

    cePtr_(nullptr),
    fePtr_(nullptr),
    pePtr_(nullptr),
    ppPtr_(nullptr),
    cpPtr_(nullptr),

    labels_(0),

    cellCentresPtr_(nullptr),
    faceCentresPtr_(nullptr),
    cellVolumesPtr_(nullptr),
    faceAreasPtr_(nullptr)
{}


Foam::primitiveMesh::primitiveMesh
(
    const label nPoints,
    const label nInternalFaces,
    const label nFaces,
    const label nCells
)
:
    nInternalPoints_(-1),
    nPoints_(nPoints),
    nEdges_(-1),
    nInternalFaces_(nInternalFaces),
    nFaces_(nFaces),
    nCells_(nCells),

    cellShapesPtr_(nullptr),
    edgesPtr_(nullptr),
    ccPtr_(nullptr),
    ecPtr_(nullptr),
    pcPtr_(nullptr),

    cfPtr_(nullptr),
    efPtr_(nullptr),
    pfPtr_(nullptr),

    cePtr_(nullptr),
    fePtr_(nullptr),
    pePtr_(nullptr),
    ppPtr_(nullptr),
    cpPtr_(nullptr),

    labels_(0),

    cellCentresPtr_(nullptr),
    faceCentresPtr_(nullptr),
    cellVolumesPtr_(nullptr),
    faceAreasPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::primitiveMesh::~primitiveMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::primitiveMesh::calcPointOrder
(
    label& nInternalPoints,
    labelList& oldToNew,
    const faceList& faces,
    const label nInternalFaces,
    const label nPoints
)
{
    // Internal points are points that are not used by a boundary face.

    // Map from old to new position
    oldToNew.setSize(nPoints);
    oldToNew = -1;


    // 1. Create compact addressing for boundary points. Start off by indexing
    // from 0 inside oldToNew. (shifted up later on)

    label nBoundaryPoints = 0;
    for (label facei = nInternalFaces; facei < faces.size(); facei++)
    {
        const face& f = faces[facei];

        forAll(f, fp)
        {
            label pointi = f[fp];

            if (oldToNew[pointi] == -1)
            {
                oldToNew[pointi] = nBoundaryPoints++;
            }
        }
    }

    // Now we know the number of boundary and internal points

    nInternalPoints = nPoints - nBoundaryPoints;

    // Move the boundary addressing up
    forAll(oldToNew, pointi)
    {
        if (oldToNew[pointi] != -1)
        {
            oldToNew[pointi] += nInternalPoints;
        }
    }


    // 2. Compact the internal points. Detect whether internal and boundary
    // points are mixed.

    label internalPointi = 0;

    bool ordered = true;

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const face& f = faces[facei];

        forAll(f, fp)
        {
            label pointi = f[fp];

            if (oldToNew[pointi] == -1)
            {
                if (pointi >= nInternalPoints)
                {
                    ordered = false;
                }
                oldToNew[pointi] = internalPointi++;
            }
        }
    }

    return ordered;
}


void Foam::primitiveMesh::reset
(
    const label nPoints,
    const label nInternalFaces,
    const label nFaces,
    const label nCells
)
{
    clearOut();

    nPoints_ = nPoints;
    nEdges_ = -1;
    nInternal0Edges_ = -1;
    nInternal1Edges_ = -1;
    nInternalEdges_ = -1;

    nInternalFaces_ = nInternalFaces;
    nFaces_ = nFaces;
    nCells_ = nCells;

    // Check if points are ordered
    label nInternalPoints;
    labelList pointMap;

    bool isOrdered = calcPointOrder
    (
        nInternalPoints,
        pointMap,
        faces(),
        nInternalFaces_,
        nPoints_
    );

    if (isOrdered)
    {
        nInternalPoints_ = nInternalPoints;
    }
    else
    {
        nInternalPoints_ = -1;
    }

    if (debug)
    {
        Pout<< "primitiveMesh::reset : mesh reset to"
            << " nInternalPoints:" << nInternalPoints_
            << " nPoints:" << nPoints_
            << " nEdges:" << nEdges_
            << " nInternalFaces:" << nInternalFaces_
            << " nFaces:" << nFaces_
            << " nCells:" << nCells_
            << endl;
    }
}


void Foam::primitiveMesh::reset
(
    const label nPoints,
    const label nInternalFaces,
    const label nFaces,
    const label nCells,
    cellList& clst
)
{
    reset
    (
        nPoints,
        nInternalFaces,
        nFaces,
        nCells
    );

    cfPtr_ = new cellList(clst, true);
}


void Foam::primitiveMesh::reset
(
    const label nPoints,
    const label nInternalFaces,
    const label nFaces,
    const label nCells,
    cellList&& clst
)
{
    reset
    (
        nPoints,
        nInternalFaces,
        nFaces,
        nCells
    );

    cfPtr_ = new cellList(move(clst));
}


Foam::tmp<Foam::scalarField> Foam::primitiveMesh::movePoints
(
    const pointField& newPoints,
    const pointField& oldPoints
)
{
    if (newPoints.size() <  nPoints() || oldPoints.size() < nPoints())
    {
        FatalErrorInFunction
            << "Cannot move points: size of given point list smaller "
            << "than the number of active points"
            << abort(FatalError);
    }

    // Create swept volumes
    const faceList& f = faces();

    tmp<scalarField> tsweptVols(new scalarField(f.size()));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(f, facei)
    {
        sweptVols[facei] = f[facei].sweptVol(oldPoints, newPoints);
    }

    // Force recalculation of all geometric data with new points
    clearGeom();

    return tsweptVols;
}


const Foam::cellShapeList& Foam::primitiveMesh::cellShapes() const
{
    if (!cellShapesPtr_)
    {
        calcCellShapes();
    }

    return *cellShapesPtr_;
}


// ************************************************************************* //
