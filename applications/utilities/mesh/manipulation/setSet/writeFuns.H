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

Class
    Foam::writeFuns

Description
    Various functions for collecting and writing binary data.

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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class writeFuns Declaration
\*---------------------------------------------------------------------------*/

class writeFuns
{
    //- Swap halves of word.
    static void swapWord(label& word32);
    static void swapWords(const label nWords, label* words32);

public:

    //- Write ascii or binary. If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, DynamicList<floatScalar>&);
    static void write(std::ostream&, const bool, DynamicList<label>&);

    //- Write ascii or binary. If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, List<floatScalar>&);
    static void write(std::ostream&, const bool, labelList&);

    //- Append point to DynamicList
    static void insert(const point&, DynamicList<floatScalar>& dest);

    //- Append elements of labelList to DynamicList
    static void insert(const labelList&, DynamicList<label>&);
    //- Append elements of scalarList to DynamicList
    static void insert(const List<scalar>&, DynamicList<floatScalar>&);

    //- Append elements of scalarList to DynamicList using map
    static void insert
    (
        const labelList& map,
        const List<scalar>& source,
        DynamicList<floatScalar>&
    );

    //- Append points to DynamicList of floats
    static void insert(const List<point>& source, DynamicList<floatScalar>&);
    //- As above but using map
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
