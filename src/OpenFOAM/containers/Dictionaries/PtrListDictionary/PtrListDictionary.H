/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::PtrListDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on the memory managed PtrList
    of \<T\> to provide ordered indexing in addition to the dictionary lookup.

SourceFiles
    PtrListDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef PtrListDictionary_H
#define PtrListDictionary_H

#include "DictionaryBase.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PtrListDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrListDictionary
:
    public DictionaryBase<PtrList<T>, T>
{

public:

    // Constructors

        //- Construct given initial list size
        PtrListDictionary(const label size);

        //- Copy construct
        PtrListDictionary(const PtrListDictionary&);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrListDictionary(Istream&, const INew&);

        //- Construct from Istream
        PtrListDictionary(Istream&);


    // Member functions

        //- Set element to pointer provided and return old element
        autoPtr<T> set(const label, const word& key, T*);

        //- Set element to autoPtr value provided and return old element
        autoPtr<T> set(const label, const word& key, autoPtr<T>&);

        //- Set element to tmp value provided and return old element
        autoPtr<T> set(const label, const word& key, tmp<T>&);


    // Member operators

        using PtrList<T>::operator[];

        //- Find and return entry
        const T& operator[](const word& key) const
        {
            return *DictionaryBase<PtrList<T>, T>::operator[](key);
        }

        //- Find and return entry
        T& operator[](const word& key)
        {
            return *DictionaryBase<PtrList<T>, T>::operator[](key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrListDictionary.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
