/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "IOstreams.H"
#include "Swap.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// return
//  -  0: different
//  - +1: identical
//  - -1: same edge, but different orientation
inline int Foam::edge::compare(const edge& a, const edge& b)
{
    if (a[0] == b[0] && a[1] == b[1])
    {
        return 1;
    }
    else if (a[0] == b[1] && a[1] == b[0])
    {
        return -1;
    }
    else
    {
        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::edge::edge()
{}


inline Foam::edge::edge(const label a, const label b)
{
    start() = a;
    end() = b;
}


inline Foam::edge::edge(const FixedList<label, 2>& a)
{
    start() = a[0];
    end() = a[1];
}


inline Foam::edge::edge(Istream& is)
:
    FixedList<label, 2>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::edge::start() const
{
    return operator[](0);
}

inline Foam::label& Foam::edge::start()
{
    return operator[](0);
}


inline Foam::label Foam::edge::end() const
{
    return operator[](1);
}

inline Foam::label& Foam::edge::end()
{
    return operator[](1);
}


inline bool Foam::edge::connected(const edge& a) const
{
    if
    (
        start() == a.start()
     || start() == a.end()
     || end() == a.start()
     || end() == a.end()
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline Foam::label Foam::edge::commonVertex(const edge& a) const
{
    if (start() == a.start() || start() == a.end())
    {
        return start();
    }
    else if (end() == a.start() || end() == a.end())
    {
        return end();
    }
    else
    {
        // No shared vertex.
        return -1;
    }
}


inline Foam::label Foam::edge::otherVertex(const label a) const
{
    if (a == start())
    {
        return end();
    }
    else if (a == end())
    {
        return start();
    }
    else
    {
        // The given vertex is not on the edge in the first place.
        return -1;
    }
}

inline void Foam::edge::flip()
{
    Swap(operator[](0), operator[](1));
}


inline Foam::edge Foam::edge::reverseEdge() const
{
    return edge(end(), start());
}


inline Foam::point Foam::edge::centre(const pointField& p) const
{
    return 0.5*(p[start()] + p[end()]);
}


inline Foam::vector Foam::edge::vec(const pointField& p) const
{
    return p[end()] - p[start()];
}


inline Foam::scalar Foam::edge::mag(const pointField& p) const
{
    return ::Foam::mag(vec(p));
}


inline Foam::linePointRef Foam::edge::line(const pointField& p) const
{
    return linePointRef(p[start()], p[end()]);
}


// * * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * //

inline bool Foam::operator==(const edge& a, const edge& b)
{
    return edge::compare(a,b) != 0;
}


inline bool Foam::operator!=(const edge& a, const edge& b)
{
    return edge::compare(a,b) == 0;
}


// ************************************************************************* //
