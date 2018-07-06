/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::faceAreaIntersect::setTriPoints
(
    const point& a,
    const point& b,
    const point& c,
    label& count,
    FixedList<triPoints, 10>& tris
) const
{
    triPoints& tp = tris[count++];
    tp[0] = a;
    tp[1] = b;
    tp[2] = c;
}


inline Foam::faceAreaIntersect::triPoints Foam::faceAreaIntersect::getTriPoints
(
    const pointField& points,
    const face& f,
    const bool reverse
) const
{
    triPoints result;

    if (reverse)
    {
        result[2] = points[f[0]];
        result[1] = points[f[1]];
        result[0] = points[f[2]];
    }
    else
    {
        result[0] = points[f[0]];
        result[1] = points[f[1]];
        result[2] = points[f[2]];
    }

    return result;
}


inline Foam::point Foam::faceAreaIntersect::planeIntersection
(
    const FixedList<scalar, 3>& d,
    const triPoints& t,
    const label negI,
    const label posI
) const
{
    return (d[posI]*t[negI] - d[negI]*t[posI])/(-d[negI] + d[posI]);
}


inline Foam::scalar Foam::faceAreaIntersect::triArea(const triPoints& t) const
{
    return mag(0.5*((t[1] - t[0])^(t[2] - t[0])));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar& Foam::faceAreaIntersect::tolerance()
{
    return tol;
}


// ************************************************************************* //
