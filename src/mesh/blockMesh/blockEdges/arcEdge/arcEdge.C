/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "arcEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(arcEdge, 0);
    addToRunTimeSelectionTable(blockEdge, arcEdge, Istream);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockEdges::arcEdge::calc(const point& pA)
{
    const tensor dps(pA - p0_, pA - p1_, (pA - p0_)^(pA - p1_));
    const tensor pMs((p0_ + pA)/2, (p1_ + pA)/2, pA);

    if (det(dps) < vSmall)
    {
        FatalErrorInFunction
            << "Arc point " << pA << " lies on the line of edge "
            << edge(start_, end_) << abort(FatalError);
    }

    const point c =
        inv(dps)
      & vector(dps.x() & pMs.x(), dps.y() & pMs.y(), dps.z() & pMs.z());

    const vector r0 = p0_ - c, r1 = p1_ - c;
    const scalar cosT = (r0 & r1)/(mag(r0)*mag(r1));
    scalar t = acos(max(-1, min(cosT, 1)));
    if (((r0 ^ r1) & dps.z()) > 0)
    {
        t = 2*constant::mathematical::pi - t;
    }

    centre_ = c;
    axis_ = - normalised(dps.z());
    theta_ = t;
    length_ = 0;
}


void Foam::blockEdges::arcEdge::calc(const scalar theta, const vector& axis)
{
    if (0 >= theta || theta >= 360)
    {
        FatalErrorInFunction
            << "Arc angle for edge " << edge(start_, end_)
            << " must take a value between 0 and 360 degrees"
            << abort(FatalError);
    }

    const vector dp = p1_ - p0_;

    const vector pM = (p0_ + p1_)/2;
    const vector rM = normalised(dp ^ axis);

    const scalar l = dp & axis;

    const vector chord = dp - l*axis;
    const scalar magChord = mag(chord);

    centre_ = pM - l*axis/2 - rM*magChord/2/tan(degToRad(theta)/2);
    axis_ = axis;
    theta_ = degToRad(theta);
    length_ = l;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::arcEdge::arcEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    p0_(points_[start_]),
    p1_(points_[end_]),
    centre_(),
    axis_(),
    theta_(),
    length_()
{
    token firstToken(is);
    is.putBack(firstToken);

    if (firstToken == token::BEGIN_LIST)
    {
        const point pA(is);
        calc(pA);
    }
    else
    {
        const scalar theta = readScalar(is);
        const vector axis(is);
        calc(theta, normalised(axis));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::arcEdge::position(const scalar lambda) const
{
    if (lambda < - small || lambda > 1 + small)
    {
        FatalErrorInFunction
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    const scalar t = theta_*lambda;
    const vector r1 = p0_ - centre_, r2 = axis_ ^ r1;

    return centre_ + r1*cos(t) + r2*sin(t) + axis_*length_*lambda;
}


Foam::scalar Foam::blockEdges::arcEdge::length() const
{
    const vector r1 = p0_ - centre_;

    // Length of a helical segment
    return degToRad(theta_*sqrt(magSqr(r1) + sqr(length_)));
}


// ************************************************************************* //
