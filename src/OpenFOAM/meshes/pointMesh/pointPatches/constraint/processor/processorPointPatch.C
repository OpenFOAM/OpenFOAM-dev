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

#include "processorPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "faceList.H"
#include "primitiveFacePatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        processorPointPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::processorPointPatch::initGeometry(PstreamBuffers& pBufs)
{
    // Algorithm:
    // Depending on whether the patch is a master or a slave, get the primitive
    // patch points and filter away the points from the global patch.

    // Create the reversed patch and pick up its points
    // so that the order is correct
    const polyPatch& pp = patch();

    faceList masterFaces(pp.size());

    forAll(pp, faceI)
    {
        masterFaces[faceI] = pp[faceI].reverseFace();
    }

    reverseMeshPoints_ = primitiveFacePatch
    (
        masterFaces,
        pp.points()
    ).meshPoints();
}


void Foam::processorPointPatch::calcGeometry(PstreamBuffers& pBufs)
{}


void Foam::processorPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::processorPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::processorPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    processorPointPatch::initGeometry(pBufs);
}


void Foam::processorPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    processorPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorPointPatch::processorPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    procPolyPatch_(refCast<const processorPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPointPatch::~processorPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::processorPointPatch::reverseMeshPoints() const
{
    return reverseMeshPoints_;
}


// ************************************************************************* //
