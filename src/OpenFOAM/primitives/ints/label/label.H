/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

Typedef
    Foam::label

Description
    A label is an int32_t or int64_t as specified by the pre-processor macro
    WM_LABEL_SIZE.

    A readLabel function is defined so that label can be constructed from
    Istream.

\*---------------------------------------------------------------------------*/

#ifndef label_H
#define label_H

#include "int.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define INT_ADD_SIZE(x,s,y) x ## s ## y
#define INT_ADD_DEF_SIZE(x,s,y) INT_ADD_SIZE(x,s,y)
#define INT_SIZE(x,y) INT_ADD_DEF_SIZE(x,WM_LABEL_SIZE,y)

#if WM_LABEL_SIZE != 32 && WM_LABEL_SIZE != 64
    #error "label.H: WM_LABEL_SIZE must be set to either 32 or 64"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef INT_SIZE(int, _t) label;

static const label labelMin = INT_SIZE(INT, _MIN);
static const label labelMax = INT_SIZE(INT, _MAX);

inline label readLabel(Istream& is)
{
    return INT_SIZE(readInt,) (is);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Raise one label to the power of another
label pow(label a, label b);

//- Return n! : 0 < n <= 12
label factorial(label n);

//- Return the log base 2 by successive bit-shifting of the given label
inline label log2(label i)
{
    label log2i = 0;
    while (i >>= 1) log2i++;
    return log2i;
}

inline label& setComponent(label& l, const direction)
{
    return l;
}

inline label component(const label l, const direction)
{
    return l;
}

inline label cmptMax(const label l)
{
    return l;
}

inline label cmptMin(const label l)
{
    return l;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelSpecific.H"

#undef INT_ADD_SIZE
#undef INT_ADD_DEF_SIZE
#undef INT_SIZE

#endif

// ************************************************************************* //
