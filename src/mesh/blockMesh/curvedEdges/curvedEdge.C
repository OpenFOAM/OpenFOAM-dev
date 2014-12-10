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
#include "curvedEdge.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(curvedEdge, 0);
    defineRunTimeSelectionTable(curvedEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvedEdge::curvedEdge
(
    const pointField& points,
    const label start,
    const label end
)
:
    points_(points),
    start_(start),
    end_(end)
{}


Foam::curvedEdge::curvedEdge(const pointField& points, Istream& is)
:
    points_(points),
    start_(readLabel(is)),
    end_(readLabel(is))
{}


Foam::curvedEdge::curvedEdge(const curvedEdge& c)
:
    points_(c.points_),
    start_(c.start_),
    end_(c.end_)
{}


Foam::autoPtr<Foam::curvedEdge> Foam::curvedEdge::clone() const
{
    notImplemented("curvedEdge::clone() const");
    return autoPtr<curvedEdge>(NULL);
}


Foam::autoPtr<Foam::curvedEdge> Foam::curvedEdge::New
(
    const pointField& points,
    Istream& is
)
{
    if (debug)
    {
        Info<< "curvedEdge::New(const pointField&, Istream&) : "
            << "constructing curvedEdge"
            << endl;
    }

    const word edgeType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(edgeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("curvedEdge::New(const pointField&, Istream&)")
            << "Unknown curvedEdge type "
            << edgeType << nl << nl
            << "Valid curvedEdge types are" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<curvedEdge>(cstrIter()(points, is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::curvedEdge::appendEndPoints
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherKnots
)
{
    pointField allKnots(otherKnots.size() + 2);

    // start/end knots
    allKnots[0] = points[start];
    allKnots[otherKnots.size() + 1] = points[end];

    // intermediate knots
    forAll(otherKnots, knotI)
    {
        allKnots[knotI+1] = otherKnots[knotI];
    }

    return allKnots;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::curvedEdge::operator=(const curvedEdge&)
{
    notImplemented("void curvedEdge::operator=(const curvedEdge&)");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const curvedEdge& p)
{
    os << p.start_ << tab << p.end_ << endl;

    return os;
}


// ************************************************************************* //
