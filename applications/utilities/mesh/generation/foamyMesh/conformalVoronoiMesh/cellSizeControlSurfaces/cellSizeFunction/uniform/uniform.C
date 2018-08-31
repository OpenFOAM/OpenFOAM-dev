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

#include "uniform.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniform, 0);
addToRunTimeSelectionTable(cellSizeFunction, uniform, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniform::uniform
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool uniform::sizeLocations
(
    const pointIndexHit& hitPt,
    const vector& n,
    pointField& shapePts,
    scalarField& shapeSizes
) const
{
    shapePts.setSize(0);
    shapeSizes.setSize(0);

    return true;
}


bool uniform::cellSize
(
    const point& pt,
    scalar& size
) const
{
    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, sqr(great)),
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

        size = 0;

        List<pointIndexHit> closeToSurfaceHits;

        surface_.findNearest
        (
            pointField(1, pt),
            scalarField(1, sqr(snapToSurfaceTol_)),
            regionIndices_,
            closeToSurfaceHits
        );

        const pointIndexHit& closeToSurface = closeToSurfaceHits[0];

        // If the nearest point is essentially on the surface, do not do a
        // getVolumeType calculation, as it will be prone to error.
        if (closeToSurface.hit())
        {
            size = surfaceCellSizeFunction_().interpolate(hitPt, index);

            return true;
        }

        pointField ptF(1, pt);
        List<volumeType> vTL(1);

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


bool uniform::setCellSize
(
    const pointField& pts
)
{
//    labelHashSet surfaceAlreadyHit(cellSize_.size());
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
//            scalarField(1, sqr(great)),
//            hits
//        );
//
//        if (hits[0].hit() && !surfaceAlreadyHit.found(hits[0].index()))
//        {
//            surfaceCellSizeFunction_().refineCellSize(hits[0].index());
//
//            surfaceAlreadyHit.insert(hits[0].index());
//        }
//    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
