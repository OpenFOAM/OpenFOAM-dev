/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "polygonTriangulate.H"

// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

inline Foam::label Foam::polygonTriangulate::renumberPolyToGlobal
(
    const label triPoly,
    const UList<label>& polyToGlobal
)
{
    return triPoly < polyToGlobal.size() ? polyToGlobal[triPoly] : -1;
}

template<class Type>
inline Type Foam::polygonTriangulate::renumberPolyToGlobal
(
    const Type& triPoly,
    const UList<label>& polyToGlobal
)
{
    Type result;
    result.resize(triPoly.size());
    forAll(triPoly, i)
    {
        result[i] = renumberPolyToGlobal(triPoly[i], polyToGlobal);
    }
    return result;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::UList<Foam::triFace>&
Foam::polygonTriangulate::triPoints() const
{
    return triPoints_;
}


inline Foam::List<Foam::triFace> Foam::polygonTriangulate::triPoints
(
    const UList<label>& polyPoints
) const
{
    return renumberPolyToGlobal(triPoints_, polyPoints);
}


inline Foam::triFace Foam::polygonTriangulate::triPoints
(
    const label trii,
    const UList<label>& polyPoints
) const
{
    return renumberPolyToGlobal(triPoints_[trii], polyPoints);
}


inline const Foam::UList<Foam::FixedList<Foam::label, 3>>&
Foam::polygonTriangulate::triEdges() const
{
    return triEdges_;
}


inline Foam::List<Foam::FixedList<Foam::label, 3>>
Foam::polygonTriangulate::triEdges
(
    const UList<label>& polyEdges
) const
{
    return renumberPolyToGlobal(triEdges_, polyEdges);
}


inline Foam::FixedList<Foam::label, 3> Foam::polygonTriangulate::triEdges
(
    const label trii,
    const UList<label>& polyEdges
) const
{
    return renumberPolyToGlobal(triEdges_[trii], polyEdges);
}


inline const Foam::UList<Foam::triFace>& Foam::polygonTriangulate::triangulate
(
    const UIndirectList<point>& points,
    const vector& normal,
    const bool simple,
    const bool optimal
)
{
    points_ = points;

    return triangulate(points_, normal, simple, optimal);
}


inline const Foam::UList<Foam::triFace>& Foam::polygonTriangulate::triangulate
(
    const UIndirectList<point>& points,
    const bool simple,
    const bool optimal
)
{
    points_ = points;

    return triangulate(points_, simple, optimal);
}


// ************************************************************************* //
