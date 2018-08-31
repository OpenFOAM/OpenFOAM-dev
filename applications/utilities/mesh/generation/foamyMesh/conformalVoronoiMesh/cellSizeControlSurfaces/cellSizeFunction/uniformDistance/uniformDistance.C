/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "uniformDistance.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniformDistance, 0);
addToRunTimeSelectionTable(cellSizeFunction, uniformDistance, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniformDistance::uniformDistance
(
    const dictionary& initialPointsDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize,
    const labelList regionIndices
)
:
    cellSizeFunction
    (
        typeName,
        initialPointsDict,
        surface,
        defaultCellSize,
        regionIndices
    ),
    distance_
    (
        readScalar(coeffsDict().lookup("distanceCoeff"))*defaultCellSize
    ),
    distanceSqr_(sqr(distance_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool uniformDistance::sizeLocations
(
    const pointIndexHit& hitPt,
    const vector& n,
    pointField& shapePts,
    scalarField& shapeSizes
) const
{
    const Foam::point& pt = hitPt.hitPoint();

    const scalar distanceCellSize =
        surfaceCellSizeFunction_().interpolate(pt, hitPt.index());

    if (sideMode_ == rmBothsides)
    {
        shapePts.resize(2);
        shapeSizes.resize(2);

        shapePts[0] = pt - n*distance_;
        shapeSizes[0] = distanceCellSize;

        shapePts[1] = pt + n*distance_;
        shapeSizes[1] = distanceCellSize;
    }
    else if (sideMode_ == smInside)
    {
        shapePts.resize(1);
        shapeSizes.resize(1);

        shapePts[0] = pt - n*distance_;
        shapeSizes[0] = distanceCellSize;
    }
    else if (sideMode_ == smOutside)
    {
        shapePts.resize(1);
        shapeSizes.resize(1);

        shapePts[0] = pt - n*distance_;
        shapeSizes[0] = distanceCellSize;
    }

    return false;
}


bool uniformDistance::cellSize
(
    const point& pt,
    scalar& size
) const
{
    size = 0;

    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, distanceSqr_),
        regionIndices_,
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    if (hitInfo.hit())
    {
        const point& hitPt = hitInfo.hitPoint();
        const label index = hitInfo.index();

        if (sideMode_ == rmBothsides)
        {
            size = surfaceCellSizeFunction_().interpolate(hitPt, index);

            return true;
        }

        // If the nearest point is essentially on the surface, do not do a
        // getVolumeType calculation, as it will be prone to error.
        if (mag(pt  - hitInfo.hitPoint()) < snapToSurfaceTol_)
        {
            size = surfaceCellSizeFunction_().interpolate(hitPt, index);

            return true;
        }

        pointField ptF(1, pt);
        List<volumeType> vTL;

        surface_.getVolumeType(ptF, vTL);

        bool functionApplied = false;

        if
        (
            sideMode_ == smInside
         && vTL[0] == volumeType::inside
        )
        {
            size = surfaceCellSizeFunction_().interpolate(hitPt, index);

            functionApplied = true;
        }
        else if
        (
            sideMode_ == smOutside
         && vTL[0] == volumeType::outside
        )
        {
            size = surfaceCellSizeFunction_().interpolate(hitPt, index);

            functionApplied = true;
        }

        return functionApplied;
    }

    return false;
}


bool uniformDistance::setCellSize
(
    const pointField& pts
)
{
//    labelHashSet surfaceAlreadyHit(surface_.size());
//
//    forAll(pts, ptI)
//    {
//        const Foam::point& pt = pts[ptI];
//
//        List<pointIndexHit> hits;
//
//        surface_.findNearest
//        (
//            pointField(1, pt),
//            scalarField(1, distanceSqr_),
//            regionIndices_,
//            hits
//        );
//
//        if (hits[0].hit() && !surfaceAlreadyHit.found(hits[0].index()))
//        {
//            surfaceCellSizeFunction_().refineSurfaceSize(hits[0].index());
//
//            surfaceAlreadyHit.insert(hits[0].index());
//        }
//    }
//
//    return true;

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
