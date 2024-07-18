/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "PtrListDictionary.H"
#include "UPtrListDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(const label size)
:
    DictionaryBase<PtrList<T>, T>(2*size)
{
    PtrList<T>::setSize(size);
}


template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(const PtrListDictionary& dict)
:
    DictionaryBase<PtrList<T>, T>(dict)
{}


template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(PtrListDictionary&& dict)
:
    DictionaryBase<PtrList<T>, T>(move(dict))
{}


template<class T>
template<class INew>
Foam::PtrListDictionary<T>::PtrListDictionary(Istream& is, const INew& iNew)
:
    DictionaryBase<PtrList<T>, T>(is, iNew)
{}


template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(Istream& is)
:
    DictionaryBase<PtrList<T>, T>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::PtrListDictionary<T>::findIndex(const word& key) const
{
    forAll(*this, i)
    {
        if (operator[](i).keyword() == key)
        {
            return i;
        }
    }

    return -1;
}


template<class T>
Foam::List<Foam::label> Foam::PtrListDictionary<T>::findIndices
(
    const wordRe& key
) const
{
    List<label> indices;

    if (!key.empty())
    {
        if (key.isPattern())
        {
            indices = findStrings(key, this->toc());
        }
        else
        {
            indices.setSize(this->size());
            label nFound = 0;
            forAll(*this, i)
            {
                if (key == operator[](i).keyword())
                {
                    indices[nFound++] = i;
                }
            }
            indices.setSize(nFound);
        }
    }

    return indices;
}


template<class T>
inline void Foam::PtrListDictionary<T>::append
(
    const word& key,
    T* ptr
)
{
    if (!DictionaryBase<PtrList<T>, T>::hashedTs_.insert(key, ptr))
    {
        FatalErrorInFunction
            << "Cannot insert with key '" << key << "' into hash-table"
            << abort(FatalError);
    }
    return PtrList<T>::append(ptr);
}


template<class T>
inline void Foam::PtrListDictionary<T>::append
(
    const word& key,
    const autoPtr<T>& aptr
)
{
    return append(key, const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline void Foam::PtrListDictionary<T>::append
(
    const word& key,
    const tmp<T>& t
)
{
    return append(key, const_cast<tmp<T>&>(t).ptr());
}


template<class T>
inline void Foam::PtrListDictionary<T>::append(T* ptr)
{
    append(ptr->keyword(), ptr);
}


template<class T>
inline void Foam::PtrListDictionary<T>::append(const autoPtr<T>& aptr)
{
    append(aptr->keyword(), aptr);
}


template<class T>
inline void Foam::PtrListDictionary<T>::append(const tmp<T>& t)
{
    append(t->keyword(), t);
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const word& key,
    T* ptr
)
{
    if (ptr == nullptr)
    {
        // Remove key from hash table if the pointer is null
        DictionaryBase<PtrList<T>, T>::hashedTs_.erase(key);
    }
    else if (!DictionaryBase<PtrList<T>, T>::hashedTs_.set(key, ptr))
    {
        FatalErrorInFunction
            << "Cannot set with key '" << key << "' into hash-table"
            << abort(FatalError);
    }

    return PtrList<T>::set(i, ptr);
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const word& key,
    const autoPtr<T>& aptr
)
{
    return set(i, key, const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const word& key,
    const tmp<T>& t
)
{
    return set(i, key, const_cast<tmp<T>&>(t).ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    T* ptr
)
{
    return set(i, ptr->keyword(), ptr);
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const autoPtr<T>& aptr
)
{
    return set(i, aptr->keyword(), aptr);
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const tmp<T>& t
)
{
    return set(i, t->keyword(), t);
}


template<class T>
template<class T2>
Foam::UPtrListDictionary<T2> Foam::PtrListDictionary<T>::convert()
{
    UPtrListDictionary<T2> result(this->size());

    forAll(*this, i)
    {
        result.set
        (
            i,
            this->operator()(i)->keyword(),
            dynamic_cast<T2*>(this->operator()(i))
        );
    }

    return result;
}


template<class T>
template<class T2>
Foam::UPtrListDictionary<T2>
Foam::PtrListDictionary<T>::lookupType()
{
    UPtrListDictionary<T2> result(this->size());

    label n = 0;
    forAll(*this, i)
    {
        if (isA<T2>(*this->operator()(i)))
        {
            result.set
            (
                n++,
                this->operator()(i)->keyword(),
                dynamic_cast<T2*>(this->operator()(i))
            );
        }
    }

    result.setSize(n);

    return result;
}


template<class T>
template<class T2>
Foam::UPtrListDictionary<const T2>
Foam::PtrListDictionary<T>::lookupType() const
{
    UPtrListDictionary<const T2> result(this->size());

    label n = 0;
    forAll(*this, i)
    {
        if (isA<T2>(*this->operator()(i)))
        {
            result.set
            (
                n++,
                this->operator()(i)->keyword(),
                dynamic_cast<const T2*>(this->operator()(i))
            );
        }
    }

    result.setSize(n);

    return result;
}


// ************************************************************************* //
