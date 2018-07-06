/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::blockEdge::start() const
{
    return start_;
}


inline Foam::label Foam::blockEdge::end() const
{
    return end_;
}


inline int Foam::blockEdge::compare(const label start, const label end) const
{
    if (start_ == start && end_ == end)
    {
        return 1;
    }
    else if (start_ == end && end_ == start)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}


inline int Foam::blockEdge::compare(const blockEdge& e) const
{
    return Foam::blockEdge::compare(e.start(), e.end());
}


inline int Foam::blockEdge::compare(const edge& e) const
{
    return Foam::blockEdge::compare(e.start(), e.end());
}


// ************************************************************************* //
