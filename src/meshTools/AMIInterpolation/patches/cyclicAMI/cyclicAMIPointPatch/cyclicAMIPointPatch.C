/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "cyclicAMIPointPatch.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicAMIPointPatch::initCalcGeometry(PstreamBuffers&)
{}


void Foam::cyclicAMIPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::cyclicAMIPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::cyclicAMIPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::cyclicAMIPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
//    cyclicAMIPointPatch::initCalcGeometry(pBufs);
}


void Foam::cyclicAMIPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
//    cyclicAMIPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicAMIPointPatch::cyclicAMIPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    cyclicAMIPolyPatch_(refCast<const cyclicAMIPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPointPatch::~cyclicAMIPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicAMIPointPatch::coupled() const
{
    return
        Pstream::parRun()
     || !this->boundaryMesh().mesh().mesh().time().processorCase();
}


// ************************************************************************* //
