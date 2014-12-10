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
#include "block.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::block::block
(
    const pointField& blockPointField,
    const curvedEdgeList& edges,
    Istream& is
)
:
    blockDescriptor(blockPointField, edges, is),
    vertices_(0),
    cells_(0),
    boundaryPatches_(0)
{}


Foam::block::block(const blockDescriptor& blockDesc)
:
    blockDescriptor(blockDesc),
    vertices_(0),
    cells_(0),
    boundaryPatches_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::block::~block()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::block::points() const
{
    if (vertices_.empty())
    {
        createPoints();
    }

    return vertices_;
}


const Foam::labelListList& Foam::block::cells() const
{
    if (cells_.empty())
    {
        createCells();
    }

    return cells_;
}


const Foam::labelListListList& Foam::block::boundaryPatches() const
{
    if (boundaryPatches_.empty())
    {
        createBoundary();
    }

    return boundaryPatches_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const block& b)
{
    os << b.points() << nl
       << b.cells() << nl
       << b.boundaryPatches() << endl;

    return os;
}


// ************************************************************************* //
