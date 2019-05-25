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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline Foam::pointField& Foam::edgeMesh::storedPoints()
{
    return points_;
}


inline Foam::edgeList& Foam::edgeMesh::storedEdges()
{
    return edges_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::edgeMesh::edgeMesh(const edgeMesh& em)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(em.points_),
    edges_(em.edges_),
    pointEdgesPtr_(nullptr)
{}


inline Foam::edgeMesh::edgeMesh(edgeMesh&& em)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(move(em.points_)),
    edges_(move(em.edges_)),
    pointEdgesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::pointField& Foam::edgeMesh::points() const
{
    return points_;
}


inline const Foam::edgeList& Foam::edgeMesh::edges() const
{
    return edges_;
}


inline const Foam::labelListList& Foam::edgeMesh::pointEdges() const
{
    if (pointEdgesPtr_.empty())
    {
        calcPointEdges();
    }
    return pointEdgesPtr_();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::edgeMesh::operator=(const edgeMesh& rhs)
{
    points_ = rhs.points_;
    edges_ = rhs.edges_;
    pointEdgesPtr_.clear();
}


void Foam::edgeMesh::operator=(edgeMesh&& rhs)
{
    points_ = move(rhs.points_);
    edges_ = move(rhs.edges_);
    pointEdgesPtr_.clear();
}


// ************************************************************************* //
