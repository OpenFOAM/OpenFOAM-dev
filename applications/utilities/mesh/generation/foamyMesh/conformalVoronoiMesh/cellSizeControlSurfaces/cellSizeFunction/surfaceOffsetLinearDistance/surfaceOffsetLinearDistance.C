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

#include "surfaceOffsetLinearDistance.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceOffsetLinearDistance, 0);
addToRunTimeSelectionTable
(
    cellSizeFunction,
    surfaceOffsetLinearDistance,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceOffsetLinearDistance::surfaceOffsetLinearDistance
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
    distanceCellSize_
    (
        readScalar(coeffsDict().lookup("distanceCellSizeCoeff"))
       *defaultCellSize
    ),
    surfaceOffset_
    (
        readScalar(coeffsDict().lookup("surfaceOffsetCoeff"))*defaultCellSize
    ),
    totalDistance_(),
    totalDistanceSqr_()
{
    if
    (
        coeffsDict().found("totalDistanceCoeff")
     || coeffsDict().found("linearDistanceCoeff")
    )
    {
        if
        (
            coeffsDict().found("totalDistanceCoeff")
         && coeffsDict().found("linearDistanceCoeff")
        )
        {
            FatalErrorInFunction
                << "totalDistanceCoeff and linearDistanceCoeff found, "
                << "specify one or other, not both."
                << nl << exit(FatalError) << endl;
        }

        if (coeffsDict().found("totalDistanceCoeff"))
        {
            totalDistance_ =
                readScalar(coeffsDict().lookup("totalDistanceCoeff"))
               *defaultCellSize;
        }
        else
        {
            totalDistance_ =
                readScalar(coeffsDict().lookup("linearDistanceCoeff"))
               *defaultCellSize
              + surfaceOffset_;
        }
    }
    else
    {
        FatalErrorInFunction
            << "totalDistanceCoeff or linearDistanceCoeff not found."
            << nl << exit(FatalError) << endl;
    }

    totalDistanceSqr_ = sqr(totalDistance_);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalar surfaceOffsetLinearDistance::sizeFunction
(
    const point& pt,
    scalar d,
    label index
) const
{
    const scalar interpolatedSize
        = surfaceCellSizeFunction_().interpolate(pt, index);

    if (d <= surfaceOffset_)
    {
        return interpolatedSize;
    }

    scalar gradient =
        (distanceCellSize_ - interpolatedSize)
       /(totalDistance_ - surfaceOffset_);

    scalar intercept = interpolatedSize - gradient*surfaceOffset_;

    return gradient*d + intercept;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool surfaceOffsetLinearDistance::sizeLocations
(
    const pointIndexHit& hitPt,
    const vector& n,
    pointField& shapePts,
    scalarField& shapeSizes
) const
{
    const Foam::point& pt = hitPt.hitPoint();

    const scalar offsetCellSize =
        surfaceCellSizeFunction_().interpolate(pt, hitPt.index());

    if (sideMode_ == rmBothsides)
    {
        shapePts.resize(4);
        shapeSizes.resize(4);

        shapePts[0] = pt - n*surfaceOffset_;
        shapeSizes[0] = offsetCellSize;
        shapePts[1] = pt - n*totalDistance_;
        shapeSizes[1] = distanceCellSize_;

        shapePts[2] = pt + n*surfaceOffset_;
        shapeSizes[2] = offsetCellSize;
        shapePts[3] = pt + n*totalDistance_;
        shapeSizes[3] = distanceCellSize_;
    }
    else if (sideMode_ == smInside)
    {
        shapePts.resize(2);
        shapeSizes.resize(2);

        shapePts[0] = pt - n*surfaceOffset_;
        shapeSizes[0] = offsetCellSize;
        shapePts[1] = pt - n*totalDistance_;
        shapeSizes[1] = distanceCellSize_;
    }
    else if (sideMode_ == smOutside)
    {
        shapePts.resize(2);
        shapeSizes.resize(2);

        shapePts[0] = pt + n*surfaceOffset_;
        shapeSizes[0] = offsetCellSize;
        shapePts[1] = pt + n*totalDistance_;
        shapeSizes[1] = distanceCellSize_;
    }

    return true;
}


bool surfaceOffsetLinearDistance::cellSize
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
        scalarField(1, totalDistanceSqr_),
        regionIndices_,
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    if (hitInfo.hit())
    {
        const point& hitPt = hitInfo.hitPoint();
        const label hitIndex = hitInfo.index();

        const scalar dist = mag(pt - hitPt);

        if (sideMode_ == rmBothsides)
        {
            size = sizeFunction(hitPt, dist, hitIndex);

            return true;
        }

        // If the nearest point is essentially on the surface, do not do a
        // getVolumeType calculation, as it will be prone to error.
        if (mag(pt  - hitInfo.hitPoint()) < snapToSurfaceTol_)
        {
            size = sizeFunction(hitPt, 0, hitIndex);

            return true;
        }

        pointField ptF(1, pt);
        List<volumeType> vTL;

        surface_.getVolumeType(ptF, vTL);

        bool functionApplied = false;

        if
        (
            sideMode_ == smInside
         && vTL[0] == volumeType::INSIDE
        )
        {
            size = sizeFunction(hitPt, dist, hitIndex);

            functionApplied = true;
        }
        else if
        (
            sideMode_ == smOutside
         && vTL[0] == volumeType::OUTSIDE
        )
        {
            size = sizeFunction(hitPt, dist, hitIndex);

            functionApplied = true;
        }

        return functionApplied;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
