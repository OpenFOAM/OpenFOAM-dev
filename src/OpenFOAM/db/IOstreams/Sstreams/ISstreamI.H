/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "ISstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ISstream::ISstream
(
    istream& is,
    const string& name,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    Istream(format, version, compression),
    name_(name),
    is_(is),
    buf_(bufInitialCapacity)
{
    if (is_.good())
    {
        setOpened();
        setGood();
    }
    else
    {
        setState(is_.rdstate());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::ISstream& Foam::ISstream::get(char& c)
{
    is_.get(c);
    setState(is_.rdstate());

    if (c == '\n')
    {
        lineNumber_++;
    }

    return *this;
}


inline int Foam::ISstream::peek()
{
    return is_.peek();
}


inline Foam::ISstream& Foam::ISstream::putback(const char& c)
{
    if (c == '\n')
    {
        lineNumber_--;
    }

    if (!is_.putback(c))
    {
        setBad();
    }

    setState(is_.rdstate());

    return *this;
}


// ************************************************************************* //
