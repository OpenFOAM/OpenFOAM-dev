/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "triIntersectLocation.H"

// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

inline Foam::label Foam::triIntersect::location::noi()
{
    return -labelMax;
}


// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

inline Foam::triIntersect::location::location
(
    const label srci,
    const label tgti
)
:
    labelPair(srci, tgti)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::triIntersect::location
Foam::triIntersect::location::srcPoint(const label srcPi)
{
    return location(srcPi, noi());
}


inline Foam::triIntersect::location
Foam::triIntersect::location::tgtPoint(const label tgtPi)
{
    return location(noi(), tgtPi);
}


inline Foam::triIntersect::location
Foam::triIntersect::location::srcTgtPoint(const label srcPi, const label tgtPi)
{
    return location(srcPi, tgtPi);
}


inline Foam::triIntersect::location
Foam::triIntersect::location::intersection(const label srcEi, const label tgtEi)
{
    return Foam::triIntersect::location(-1 - srcEi, -1 - tgtEi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::triIntersect::location::location()
:
    labelPair(noi(), noi())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::triIntersect::location::isNull() const
{
    return first() == noi() && second() == noi();
}


inline bool Foam::triIntersect::location::isSrcPoint() const
{
    return first() != noi() && first() >= 0;
}


inline bool Foam::triIntersect::location::isTgtPoint() const
{
    return second() != noi() && second() >= 0;
}


inline bool Foam::triIntersect::location::isSrcNotTgtPoint() const
{
    return isSrcPoint() && !isTgtPoint();
}


inline bool Foam::triIntersect::location::isTgtNotSrcPoint() const
{
    return !isSrcPoint() && isTgtPoint();
}


inline bool Foam::triIntersect::location::isSrcAndTgtPoint() const
{
    return isSrcPoint() && isTgtPoint();
}


inline bool Foam::triIntersect::location::isIntersection() const
{
    return
        first() != noi() && first() < 0
     && second() != noi() && second() < 0;
}


inline Foam::label Foam::triIntersect::location::srcPointi() const
{
    #ifdef FULLDEBUG
    if (!isSrcPoint())
    {
        FatalErrorInFunction
            << "Source point requested from non-source-point location"
            << exit(FatalError);
    }
    #endif
    return first();
}


inline Foam::label Foam::triIntersect::location::tgtPointi() const
{
    #ifdef FULLDEBUG
    if (!isTgtPoint())
    {
        FatalErrorInFunction
            << "Target point requested from non-source-point location"
            << exit(FatalError);
    }
    #endif
    return second();
}


inline Foam::label Foam::triIntersect::location::srcEdgei() const
{
    #ifdef FULLDEBUG
    if (!isIntersection())
    {
        FatalErrorInFunction
            << "Source edge requested from non-intersection location"
            << exit(FatalError);
    }
    #endif
    return -1 - first();
}


inline Foam::label Foam::triIntersect::location::tgtEdgei() const
{
    #ifdef FULLDEBUG
    if (!isIntersection())
    {
        FatalErrorInFunction
            << "Target edge requested from non-intersection location"
            << exit(FatalError);
    }
    #endif
    return -1 - second();
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

inline bool Foam::triIntersect::operator==(const location& a, const location& b)
{
    return
        static_cast<const labelPair&>(a)
     == static_cast<const labelPair&>(b);
}


inline bool Foam::triIntersect::operator!=(const location& a, const location& b)
{
    return !(a == b);
}


// ************************************************************************* //
