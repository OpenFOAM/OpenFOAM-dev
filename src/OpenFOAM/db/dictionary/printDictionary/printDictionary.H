/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::printDictionary

Description
    Enables the printing of a dictionary
    and subsequently looked-up defaulted entries

SourceFiles
    printDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef printDictionary_H
#define printDictionary_H

#include "HashSet.H"
#include "tmpNrc.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                          Class printDictionary Declaration
\*---------------------------------------------------------------------------*/

class printDictionary
{
private:

    // Private Static Member Functions

        //- Remove the default dictionary for the given dictionary
        static void removeDefaults
        (
            const dictionary* dictPtr,
            HashSet<const dictionary*, Hash<void*>>& removeDictPtrs
        );

        //- Remove the default dictionary for the given dictionary
        static void removeDefaults(const dictionary* dictPtr);

        //- Allocate the default dictionary for the given dictionary
        static void setDefaults(const dictionary& dict);

        //- Set the default dictionary pointers for the sub-dictionaries
        //  of the given dictionary
        static void setSubDefaults
        (
            const dictionary& dict,
            dictionary& defaults
        );

        //- Print the given dictionary and its accumulated defaults
        static void print(const dictionary& dict, const dictionary& defaults);


    // Private Static Data

        //- Table from the names of dictionaries to be printed to pointers
        //  to those dictionaries and the number of times they are currently
        //  marked to be printed
        static HashTable
        <
            Tuple2<const dictionary*, label>,
            fileName,
            Hash<fileName>
        > dictNameToDictPtrAndCount_;

        //- Table from pointers to the dictionaries to be printed to
        //  dictionaries containing all entries including defaults
        static HashTable
        <
            tmpNrc<dictionary>,
            const dictionary*,
            Hash<void*>
        > dictPtrToDefaults_;

        //- The dictionaries to be printed
        UPtrList<const dictionary> dicts_;

        //- Names of the dictionaries to be printed
        List<fileName> dictNames_;


    // Private Constructors

        //- Construct null
        printDictionary();


    // Private Member Functions

        //- Add a dictionary to the list to be printed
        void add(const dictionary& dict);

        //- Add a dictionary name to the list to be printed
        void add(const fileName& dictName);


public:

    // Constructors

        //- Construct adding a dictionary to the list to be printed
        template<class ... Args>
        printDictionary(const dictionary& dict, const Args& ... args)
        :
            printDictionary(args ...)
        {
            add(dict);
        }

        //- Construct adding a dictionary name to the list to be printed
        template<class ... Args>
        printDictionary(const fileName& dictName, const Args& ... args)
        :
            printDictionary(args ...)
        {
            add(dictName);
        }


    //- Destructor
    ~printDictionary();


    // Member Functions

        //- Set the dictionary to be printed
        static void set(const dictionary& dict);

        //- Return if a dictionary is marked to be printed
        static bool prints(const dictionary& dict)
        {
            return dictNameToDictPtrAndCount_.found(dict.name());
        }

        //- Return if a dictionary is marked to be printed
        static bool prints(const fileName& dictName)
        {
            return dictNameToDictPtrAndCount_.found(dictName);
        }

        //- Return if a dictionary exists to add defaults to for a given
        //  dictionary
        static bool haveDefaults(const dictionary& dict)
        {
            return dictPtrToDefaults_.found(&dict);
        }

        //- Return the dictionary to add defaults to for a given dictionary
        static dictionary& defaults(const dictionary& dict)
        {
            return const_cast<dictionary&>(dictPtrToDefaults_[&dict]());
        }
};


template<class Type, class... Args>
inline Type PrintDictionary(const dictionary& dict, Args&... args)
{
    printDictionary print(dict);
    return Type(args...);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
