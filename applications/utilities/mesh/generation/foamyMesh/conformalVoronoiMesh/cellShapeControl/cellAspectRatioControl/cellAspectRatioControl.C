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

#include "cellAspectRatioControl.H"
#include "vectorTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellAspectRatioControl::cellAspectRatioControl
(
    const dictionary& motionDict
)
:
    aspectRatioDict_(motionDict.subOrEmptyDict("cellAspectRatioControl")),
    aspectRatio_(aspectRatioDict_.lookupOrDefault<scalar>("aspectRatio", 1.0)),
    aspectRatioDirection_
    (
        aspectRatioDict_.lookupOrDefault<vector>
        (
            "aspectRatioDirection",
            vector::zero
        )
    )
{
    // Normalise the direction
    aspectRatioDirection_ /= mag(aspectRatioDirection_) + SMALL;

    Info<< nl
        << "Cell Aspect Ratio Control" << nl
        << "    Ratio     : " << aspectRatio_ << nl
        << "    Direction : " << aspectRatioDirection_
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellAspectRatioControl::~cellAspectRatioControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellAspectRatioControl::updateCellSizeAndFaceArea
(
    vector& alignmentDir,
    scalar& targetFaceArea,
    scalar& targetCellSize
) const
{
    const scalar cosAngle =
        mag(vectorTools::cosPhi(alignmentDir, aspectRatioDirection_));

    // Change target face area based on aspect ratio
    targetFaceArea +=
        targetFaceArea
       *(aspectRatio_ - 1.0)
       *(1.0 - cosAngle);

    // Change target cell size based on aspect ratio
    targetCellSize +=
        targetCellSize
       *(aspectRatio_ - 1.0)
       *cosAngle;

    alignmentDir *= 0.5*targetCellSize;
}


void Foam::cellAspectRatioControl::updateDeltaVector
(
    const vector& alignmentDir,
    const scalar targetCellSize,
    const scalar rABMag,
    vector& delta
) const
{
    const scalar cosAngle =
        mag(vectorTools::cosPhi(alignmentDir, aspectRatioDirection_));

    delta +=
        0.5
       *delta
       *cosAngle
       *(targetCellSize/rABMag)
       *(aspectRatio_ - 1.0);
}


// ************************************************************************* //
