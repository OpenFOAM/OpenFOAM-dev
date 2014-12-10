/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "error.H"
#include "lineEdge.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lineEdge, 0);
    addToRunTimeSelectionTable(curvedEdge, lineEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineEdge::lineEdge
(
    const pointField& points,
    const label start,
    const label end
)
:
    curvedEdge(points, start, end)
{}


Foam::lineEdge::lineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is)
{}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::lineEdge::~lineEdge()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::lineEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("lineEdge::position(const scalar)")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    return points_[start_] + lambda * (points_[end_] - points_[start_]);
}


Foam::scalar Foam::lineEdge::length() const
{
    return mag(points_[end_] - points_[start_]);
}


// ************************************************************************* //
