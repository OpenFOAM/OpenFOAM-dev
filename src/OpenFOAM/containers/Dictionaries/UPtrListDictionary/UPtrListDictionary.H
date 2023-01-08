/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::UPtrListDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on the UPtrList of \<T\> to
    provide ordered indexing in addition to the dictionary lookup.

SourceFiles
    UPtrListDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef UPtrListDictionary_H
#define UPtrListDictionary_H

#include "DictionaryBase.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
class PtrListDictionary;

/*---------------------------------------------------------------------------*\
                      Class UPtrListDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrListDictionary
:
    public DictionaryBase<UPtrList<T>, T>
{

public:

    // Constructors

        //- Construct given initial list size
        UPtrListDictionary(const label size);

        //- Copy constructor
        UPtrListDictionary(const UPtrListDictionary&);

        //- Copy constructor of PtrListDictionary
        UPtrListDictionary(const PtrListDictionary<T>&);

        //- Move constructor
        UPtrListDictionary(UPtrListDictionary&&);


    // Member Functions

        //- Is element set
        inline bool set(const label) const;

        //- Set element. Return old element (can be nullptr).
        //  No checks on new element
        inline T* set(const label, const word& key, T*);


    // Member Operators

        using UPtrList<T>::operator[];

        //- Find and return entry
        const T& operator[](const word& key) const
        {
            return *DictionaryBase<UPtrList<T>, T>::operator[](key);
        }

        //- Find and return entry
        T& operator[](const word& key)
        {
            return *DictionaryBase<UPtrList<T>, T>::operator[](key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPtrListDictionary.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
