/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "linearSpatial.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linearSpatial, 0);
addToRunTimeSelectionTable(cellSizeFunction, linearSpatial, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearSpatial::linearSpatial
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
    referencePoint_(coeffsDict().lookup("referencePoint")),
    referenceCellSize_
    (
        readScalar(coeffsDict().lookup("referenceCellSizeCoeff"))
       *defaultCellSize
    ),
    direction_(coeffsDict().lookup("direction")),
    cellSizeGradient_(readScalar(coeffsDict().lookup("cellSizeGradient")))
{
    direction_ /= mag(direction_);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalar linearSpatial::sizeFunction(const point& pt) const
{
    return
        referenceCellSize_
      + ((pt - referencePoint_) & direction_)*cellSizeGradient_;
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool linearSpatial::sizeLocations
(
    const pointIndexHit& hitPt,
    const vector& n,
    pointField& shapePts,
    scalarField& shapeSizes
) const
{
    if (sideMode_ == rmBothsides)
    {
    }
    else if (sideMode_ == smInside)
    {
    }
    else if (sideMode_ == smOutside)
    {
    }

    return false;
}


bool linearSpatial::cellSize
(
    const point& pt,
    scalar& size
) const
{
    if (sideMode_ == rmBothsides)
    {
        size = sizeFunction(pt);

        return true;
    }

    size = 0;

    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, sqr(snapToSurfaceTol_)),
        regionIndices_,
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    // If the nearest point is essentially on the surface, do not do a
    // getVolumeType calculation, as it will be prone to error.
    if (hitInfo.hit())
    {
        size = sizeFunction(pt);

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
        size = sizeFunction(pt);

        functionApplied = true;
    }
    else if
    (
        sideMode_ == smOutside
     && vTL[0] == volumeType::OUTSIDE
    )
    {
        size = sizeFunction(pt);

        functionApplied = true;
    }

    return functionApplied;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
