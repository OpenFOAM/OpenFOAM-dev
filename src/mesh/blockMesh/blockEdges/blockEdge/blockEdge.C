/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "blockEdge.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockEdge, 0);
    defineRunTimeSelectionTable(blockEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdge::blockEdge
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


Foam::blockEdge::blockEdge
(
    const pointField& points,
    Istream& is
)
:
    points_(points),
    start_(readLabel(is)),
    end_(readLabel(is))
{}


Foam::autoPtr<Foam::blockEdge> Foam::blockEdge::clone() const
{
    NotImplemented;
    return autoPtr<blockEdge>(nullptr);
}


Foam::autoPtr<Foam::blockEdge> Foam::blockEdge::New
(
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction << "Constructing blockEdge" << endl;
    }

    const word edgeType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(edgeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown blockEdge type "
            << edgeType << nl << nl
            << "Valid blockEdge types are" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<blockEdge>(cstrIter()(geometry, points, is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::blockEdge::appendEndPoints
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherKnots
)
{
    pointField allKnots(otherKnots.size() + 2);

    // Start/end knots
    allKnots[0] = points[start];
    allKnots[otherKnots.size() + 1] = points[end];

    // Intermediate knots
    forAll(otherKnots, knotI)
    {
        allKnots[knotI+1] = otherKnots[knotI];
    }

    return allKnots;
}


Foam::tmp<Foam::pointField>
Foam::blockEdge::position(const scalarList& lambdas) const
{
    tmp<pointField> tpoints(new pointField(lambdas.size()));
    pointField& points = tpoints.ref();

    forAll(lambdas, i)
    {
        points[i] = position(lambdas[i]);
    }
    return tpoints;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockEdge& p)
{
    os << p.start_ << tab << p.end_ << endl;

    return os;
}


// ************************************************************************* //
