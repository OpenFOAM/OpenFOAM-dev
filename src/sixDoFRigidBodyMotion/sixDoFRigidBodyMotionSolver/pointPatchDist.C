/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "pointPatchDist.H"
#include "externalPointEdgePoint.H"
#include "pointMesh.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointPatchDist::pointPatchDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const pointField& points
)
:
    pointScalarField
    (
        IOobject
        (
            "pointDistance",
            pMesh.db().time().timeName(),
            pMesh.db()
        ),
        pMesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    points_(points),
    patchIDs_(patchIDs),
    nUnset_(0)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointPatchDist::~pointPatchDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointPatchDist::correct()
{
    const pointBoundaryMesh& pbm = mesh().boundary();

    label nPoints = 0;

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchI = iter.key();
        nPoints += pbm[patchI].meshPoints().size();
    }

    externalPointEdgePoint::trackingData td(points_);

    // Set initial changed points to all the patch points(if patch present)
    List<externalPointEdgePoint> wallInfo(nPoints);
    labelList wallPoints(nPoints);
    nPoints = 0;

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchI = iter.key();
        // Retrieve the patch now we have its index in patches.

        const labelList& mp = pbm[patchI].meshPoints();

        forAll(mp, ppI)
        {
            label meshPointI = mp[ppI];
            wallPoints[nPoints] = meshPointI;
            wallInfo[nPoints] = externalPointEdgePoint
            (
                td.points_[meshPointI],
                0.0
            );
            nPoints++;
        }
    }

    // Current info on points
    List<externalPointEdgePoint> allPointInfo(mesh()().nPoints());

    // Current info on edges
    List<externalPointEdgePoint> allEdgeInfo(mesh()().nEdges());

    PointEdgeWave
    <
        externalPointEdgePoint,
        externalPointEdgePoint::trackingData
    > wallCalc
    (
        mesh()(),
        wallPoints,
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints(), // max iterations
        td
    );

    pointScalarField& psf = *this;


    forAll(allPointInfo, pointI)
    {
        if (allPointInfo[pointI].valid(td))
        {
            psf[pointI] = Foam::sqrt(allPointInfo[pointI].distSqr());
        }
        else
        {
            nUnset_++;
        }
    }
}


// ************************************************************************* //
