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

#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline const Type& operator+(const Type& t, const zero&)
{
    return t;
}

template<class Type>
inline const Type& operator+(const zero&, const Type& t)
{
    return t;
}

template<class Type>
inline const Type& operator-(const Type& t, const zero&)
{
    return t;
}

template<class Type>
inline Type operator-(const zero&, const Type& t)
{
    return -t;
}

template<class Type>
inline zero operator*(const Type& t, const zero&)
{
    return Zero;
}

template<class Type>
inline zero operator*(const zero&, const Type& t)
{
    return Zero;
}

template<class Type>
inline zero operator/(const zero&, const Type& t)
{
    return Zero;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
