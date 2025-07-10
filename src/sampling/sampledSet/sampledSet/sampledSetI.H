/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "sampledSet.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::sampledSet::name() const
{
    return name_;
}


inline const Foam::polyMesh& Foam::sampledSet::mesh() const
{
    return mesh_;
}


inline const Foam::coordSet& Foam::sampledSet::coords() const
{
    if (!coordsPtr_.valid())
    {
        setSamples();
    }
    return coordsPtr_();
}


inline Foam::label Foam::sampledSet::size() const
{
    return coords().size();
}


inline const Foam::labelList& Foam::sampledSet::segments() const
{
    return coords().segments();
}


inline const Foam::pointField& Foam::sampledSet::positions() const
{
    return coords().positions();
}


inline const Foam::labelList& Foam::sampledSet::cells() const
{
    if (!cellsPtr_.valid())
    {
        setSamples();
    }
    return cellsPtr_;
}


inline const Foam::labelList& Foam::sampledSet::faces() const
{
    if (!facesPtr_.valid())
    {
        setSamples();
    }
    return facesPtr_();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline Foam::sampledSet::operator const coordSet&() const
{
    return coords();
}


// ************************************************************************* //
