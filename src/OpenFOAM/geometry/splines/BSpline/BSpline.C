/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "BSpline.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::BSpline::position
(
    const label segmenti,
    const scalar segmentLambdaArg
) const
{
    // Quick return if segment index is out of range
    if (segmenti < 0) return start();
    if (segmenti >= nSegments()) return end();

    // Points bounding the segment
    const point& p0 = knots()[segmenti];
    const point& p1 = knots()[segmenti + 1];

    // Points before and after the segment. If at the ends, use reflection.
    const point pMinus1 =
        segmenti == 0 ? 2*p0 - p1 : knots()[segmenti - 1];
    const point p2 =
        segmenti == nSegments() - 1 ? 2*p1 - p0 : knots()[segmenti + 2];

    // Clip the argument
    const scalar segmentLambda =
        min(max(segmentLambdaArg, scalar(0)), scalar(1));

    // Evaluate the function and return
    const vector A = -pMinus1 + 3*p0 - 3*p1 + p2;
    const vector B = 3*pMinus1 - 6*p0 + 3*p1;
    const vector C = -3*pMinus1 + 3*p1;
    const vector D = pMinus1 + 4*p0 + p1;
    return (segmentLambda*(segmentLambda*(segmentLambda*A + B) + C) + D)/6;
}


// ************************************************************************* //
