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

#include "writeFuns.H"

#if defined(__mips) && !defined(__SICORTEX__)
#include <standards.h>
#include <sys/endian.h>
#endif

#if defined(LITTLE_ENDIAN) \
 || defined(_LITTLE_ENDIAN) \
 || defined(__LITTLE_ENDIAN)
#   define LITTLEENDIAN 1
#elif defined(BIG_ENDIAN) || defined(_BIG_ENDIAN) || defined(__BIG_ENDIAN)
#   undef LITTLEENDIAN
#else
#   error "Cannot find LITTLE_ENDIAN or BIG_ENDIAN symbol defined."
#   error "Please add to compilation options"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::writeFuns::swapWord(label& word32)
{
    char* mem =  reinterpret_cast<char*>(&word32);

    char a = mem[0];
    mem[0] = mem[3];
    mem[3] = a;

    a = mem[1];
    mem[1] = mem[2];
    mem[2] = a;
}


void Foam::writeFuns::swapWords(const label nWords, label* words32)
{
    for (label i = 0; i < nWords; i++)
    {
        swapWord(words32[i]);
    }
}


void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    List<floatScalar>& fField
)
{
    if (binary)
    {
#       ifdef LITTLEENDIAN
        swapWords(fField.size(),  reinterpret_cast<label*>(fField.begin()));
#       endif

        os.write
        (
            reinterpret_cast<char*>(fField.begin()),
            fField.size()*sizeof(float)
        );

        os  << std::endl;
    }
    else
    {
        forAll(fField, i)
        {
            os  << fField[i] << ' ';

            if (i > 0 && (i % 10) == 0)
            {
                os  << std::endl;
            }
        }
        os  << std::endl;
    }
}


void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    DynamicList<floatScalar>& fField
)
{
    List<floatScalar>& fld = fField.shrink();

    write(os, binary, fld);
}


void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    labelList& elems
)
{
    if (binary)
    {
#       ifdef LITTLEENDIAN
        swapWords(elems.size(),  reinterpret_cast<label*>(elems.begin()));
#       endif
        os.write
        (
            reinterpret_cast<char*>(elems.begin()),
            elems.size()*sizeof(label)
        );

        os  << std::endl;
    }
    else
    {
        forAll(elems, i)
        {
            os  << elems[i] << ' ';

            if (i > 0 && (i % 10) == 0)
            {
                os  << std::endl;
            }
        }
        os  << std::endl;
    }
}


void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    DynamicList<label>& elems
)
{
    labelList& fld = elems.shrink();

    write(os, binary, fld);
}


// Store vector in dest.
void Foam::writeFuns::insert(const point& pt, DynamicList<floatScalar>& dest)
{
    dest.append(float(pt.x()));
    dest.append(float(pt.y()));
    dest.append(float(pt.z()));
}


// Store labelList in dest.
void Foam::writeFuns::insert(const labelList& source, DynamicList<label>& dest)
{
    forAll(source, i)
    {
        dest.append(source[i]);
    }
}


// Store scalarField in dest
void Foam::writeFuns::insert
(
    const List<scalar>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
        dest.append(float(source[i]));
    }
}


// Store scalarField (indexed through map) in dest
void Foam::writeFuns::insert
(
    const labelList& map,
    const List<scalar>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(map, i)
    {
        dest.append(float(source[map[i]]));
    }
}


void Foam::writeFuns::insert
(
    const List<point>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
       insert(source[i], dest);
    }
}

void Foam::writeFuns::insert
(
    const labelList& map,
    const List<point>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(map, i)
    {
       insert(source[map[i]], dest);
    }
}


// ************************************************************************* //
