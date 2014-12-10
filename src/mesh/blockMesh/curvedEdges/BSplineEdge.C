/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "BSplineEdge.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BSplineEdge, 0);

    addToRunTimeSelectionTable
    (
        curvedEdge,
        BSplineEdge,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BSplineEdge::BSplineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& internalPoints
)
:
    curvedEdge(points, start, end),
    BSpline(appendEndPoints(points, start, end, internalPoints))
{}


Foam::BSplineEdge::BSplineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    BSpline(appendEndPoints(points, start_, end_, pointField(is)))
{
    token t(is);
    is.putBack(t);

    // discard unused start/end tangents
    if (t == token::BEGIN_LIST)
    {
        vector tangent0Ignored(is);
        vector tangent1Ignored(is);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BSplineEdge::~BSplineEdge()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::BSplineEdge::position(const scalar mu) const
{
    return BSpline::position(mu);
}


Foam::scalar Foam::BSplineEdge::length() const
{
    return BSpline::length();
}


// ************************************************************************* //
