/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::hashedWordList

Description
    A wordList with hashed indices for faster lookup by name.

SourceFiles
    hashedWordListI.H
    hashedWordList.C

\*---------------------------------------------------------------------------*/

#ifndef hashedWordList_H
#define hashedWordList_H

#include "wordList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class hashedWordList;

// Forward declaration of friend functions and operators
Istream& operator>>(Istream&, hashedWordList&);
Ostream& operator<<(Ostream&, const hashedWordList&);


/*---------------------------------------------------------------------------*\
                       Class hashedWordList Declaration
\*---------------------------------------------------------------------------*/

class hashedWordList
:
    public List<word>
{
    // Private Data

        HashTable<label, word> indices_;


    // Private Member Functions

        //- Rebuild the hash of indices
        void rehash();

public:

    // Constructors

        //- Construct null
        hashedWordList();

        //- Copy constructor
        hashedWordList(const hashedWordList&);

        //- Move constructor
        hashedWordList(hashedWordList&&);

        //- Construct from list of names
        hashedWordList(const UList<word>& names);

        //- Move constructor transferring the parameter contents
        hashedWordList(List<word>&& names);

        //- Construct from number and list of names
        hashedWordList(const label nNames, const char** names);

        //- Construct from a nullptr-terminated list of names
        hashedWordList(const char** names);

        //- Construct from Istream
        hashedWordList(Istream&);


    // Member Functions

        //- Clear the list, i.e. set size to zero.
        void clear();

        //- Append an element at the end of the list
        void append(const word&);

        //- Does the list contain the specified name
        inline bool found(const word&) const;

        //- Transfer the contents of the argument List into this list
        //  and annul the argument list.
        void transfer(List<word>&);


    // Member Operators

        //- Assignment operator
        inline void operator=(const hashedWordList&);

        //- Move assignment operator
        inline void operator=(hashedWordList&&);

        //- Assignment operator from list of names
        inline void operator=(const UList<word>& names);

        //- Return name corresponding to specified index
        inline const word& operator[](const label index) const;

        //- Return index corresponding to specified name
        inline label operator[](const word&) const;


    // Istream operators

        friend Istream& operator>>(Istream&, hashedWordList&);
        friend Ostream& operator<<(Ostream&, const hashedWordList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hashedWordListI.H"

#endif

// ************************************************************************* //
