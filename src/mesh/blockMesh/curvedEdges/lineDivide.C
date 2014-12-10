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

#include "lineDivide.H"
#include "curvedEdge.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineDivide::lineDivide
(
    const curvedEdge& cedge,
    const label ndiv,
    const scalar xratio
)
:
    points_(ndiv + 1),
    divisions_(ndiv + 1)
{
    divisions_[0]    = 0.0;
    divisions_[ndiv] = 1.0;

    // calculate the spacing
    if (xratio == 1.0)
    {
        for (label i=1; i < ndiv; i++)
        {
            divisions_[i] = scalar(i)/ndiv;
        }
    }
    else
    {
        for (label i=1; i < ndiv; i++)
        {
            divisions_[i] = (1.0 - pow(xratio, i))/(1.0 - pow(xratio, ndiv));
        }
    }

    // calculate the points
    for (label i=0; i <= ndiv; i++)
    {
        points_[i] = cedge.position(divisions_[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::lineDivide::points() const
{
    return points_;
}


const Foam::scalarList& Foam::lineDivide::lambdaDivisions() const
{
    return divisions_;
}


// ************************************************************************* //
