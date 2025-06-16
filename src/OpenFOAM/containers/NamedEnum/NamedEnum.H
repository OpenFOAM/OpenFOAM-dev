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

Class
    Foam::NamedEnum

Description
    Initialise the NamedEnum HashTable from the static list of names.

SourceFiles
    NamedEnum.C

\*---------------------------------------------------------------------------*/

#ifndef NamedEnum_H
#define NamedEnum_H

#include "FixedList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                          Class NamedEnum Declaration
\*---------------------------------------------------------------------------*/

template<class Enum, unsigned int nEnum>
class NamedEnum
:
    public FixedList<word, nEnum>
{
    // Private Member Data

        //- Map of name -> enumeration
        HashTable<unsigned int> table_;


public:

    typedef FixedList<word, nEnum> namesType;

    // Constructors

        //- Construct from list of names
        NamedEnum(std::initializer_list<word>);

        //- Disallow default bitwise copy construction
        NamedEnum(const NamedEnum&) = delete;


    // Member Functions

        //- Read a word from Istream and return the corresponding
        //  enumeration element
        Enum read(Istream&) const;

        //- Lookup the name from the dictionary if present and return the
        //  corresponding enumeration element, otherwise return the defaultValue
        Enum lookupOrDefault
        (
            const word& name,
            const dictionary& dict,
            const Enum defaultValue
        ) const;

        //- Select the name from the dictionary and return the
        //  corresponding enumeration element.
        //  If more than one of the options is present in the dictionary
        //  generate a FatalError
        Enum select(const dictionary& dict) const;

        //- Write the name representation of the enumeration to an Ostream
        void write(const Enum e, Ostream&) const;


    // Member Operators

        //- Return the enumeration element corresponding to the given name
        const Enum operator[](const char* name) const
        {
            return Enum(table_[name]);
        }

        //- Return the enumeration element corresponding to the given name
        const Enum operator[](const word& name) const
        {
            return Enum(table_[name]);
        }

        //- Return the name of the given enumeration element
        const word& operator[](const Enum e) const
        {
            return FixedList<word, nEnum>::operator[](unsigned(e));
        }

        //- Disallow default bitwise assignment
        void operator=(const NamedEnum&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NamedEnum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
