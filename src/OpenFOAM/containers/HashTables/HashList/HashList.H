/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::HashList

Description
    HashList class. Like HashTable, but much less dynamic memory-y. Should be
    faster for small sets of non-dynamic primitive types (labels, edges,
    points, etc...). It is also much less functional at present. There is no
    re-sizing, so you have to make sure it is constructed sufficiently large to
    hold all the data that will ever be inserted into it.

SourceFiles
    HashList.C

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "Tuple2.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef HashList_H
#define HashList_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HashList Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Key=word, class Hash=string::hash>
class HashList
:
    private List<Tuple2<Key, Type>>
{
public:

    // Public Static Member Data

        //- Null key value for unset elements in the list
        static const Key nullKey;


    // Constructors

        //- Construct given a size
        HashList(const label size);


    // Member Functions

        //- The size of the underlying list
        inline label capacity() const;

        //- Clear all elements
        void clear();

        //- Resize and clear all elements
        void resizeAndClear(const label newSize);

        //- Insert into the hash list. Return true if the value was newly
        //  inserted, or false if it was already there.
        bool insert(const Key& k, const Type& t);


    // Member Operators

        //- Retrieve from the hash list
        const Type& operator[](const Key& k) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HashList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
