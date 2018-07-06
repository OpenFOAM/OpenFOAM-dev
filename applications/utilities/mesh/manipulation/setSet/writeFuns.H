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

Class
    Foam::writeFuns

Description
    Various functions for collecting and writing binary data.

    The LITTLE_ENDIAN is based on 32bit words.
    It is not clear how 64bit labels should be handled, currently they are
    split into two 32bit words and swapWord applied to these two.

    writeFuns should be a namespace rather than a class.

SourceFiles
    writeFuns.C

\*---------------------------------------------------------------------------*/

#ifndef writeFuns_H
#define writeFuns_H

#include "labelList.H"
#include "floatScalar.H"
#include "OFstream.H"
#include "DynamicList.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class writeFuns Declaration
\*---------------------------------------------------------------------------*/

class writeFuns
{
    // Private member functions

        //- Swap halves of word
        static void swapWord(int32_t& word32);

        //- Swap halves of word
        static void swapWords(const label nWords, int32_t* words32);


public:

    //- Write floats ascii or binary.
    //  If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, DynamicList<floatScalar>&);

    //- Write labels ascii or binary.
    //  If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, DynamicList<label>&);

    //- Write floats ascii or binary.
    //  If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, List<floatScalar>&);

    //- Write labels ascii or binary.
    //  If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, labelList&);

    //- Append point to given DynamicList
    static void insert(const point&, DynamicList<floatScalar>& dest);

    //- Append elements of labelList to given DynamicList
    static void insert(const labelList&, DynamicList<label>&);

    //- Append elements of scalarList to given DynamicList
    static void insert(const List<scalar>&, DynamicList<floatScalar>&);

    //- Append elements of scalarList to given DynamicList using map
    static void insert
    (
        const labelList& map,
        const List<scalar>& source,
        DynamicList<floatScalar>&
    );

    //- Append points to given DynamicList of floats
    static void insert(const List<point>& source, DynamicList<floatScalar>&);

    //- Append points to given DynamicList of floats using map
    static void insert
    (
        const labelList& map,
        const List<point>& source,
        DynamicList<floatScalar>&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
