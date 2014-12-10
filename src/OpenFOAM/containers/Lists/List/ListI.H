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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T>::List()
{}


template<class T>
inline Foam::autoPtr<Foam::List<T> > Foam::List<T>::clone() const
{
    return autoPtr<List<T> >(new List<T>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::List<T>& Foam::List<T>::null()
{
    return *reinterpret_cast< List<T>* >(0);
}


template<class T>
inline void Foam::List<T>::resize(const label newSize)
{
    this->setSize(newSize);
}


template<class T>
inline void Foam::List<T>::resize(const label newSize, const T& a)
{
    this->setSize(newSize, a);
}


template<class T>
inline T& Foam::List<T>::newElmt(const label i)
{
    if (i >= this->size())
    {
        setSize(2*this->size());
    }

    return UList<T>::operator[](i);
}


template<class T>
inline void Foam::List<T>::size(const label n)
{
    UList<T>::size_ = n;
}


template<class T>
inline Foam::label Foam::List<T>::size() const
{
    return UList<T>::size_;
}


template<class T>
inline Foam::Xfer<Foam::List<T> > Foam::List<T>::xfer()
{
    return xferMove(*this);
}


template<class T>
inline void Foam::List<T>::append(const T& t)
{
    setSize(size()+1, t);
}


template<class T>
inline void Foam::List<T>::append(const UList<T>& lst)
{
    if (this == &lst)
    {
        FatalErrorIn
        (
            "List<T>::append(const UList<T>&)"
        )   << "attempted appending to self" << abort(FatalError);
    }

    label nextFree = this->size();
    setSize(nextFree + lst.size());

    forAll(lst, elemI)
    {
        this->operator[](nextFree++) = lst[elemI];
    }
}


template<class T>
inline void Foam::List<T>::append(const UIndirectList<T>& lst)
{
    label nextFree = this->size();
    setSize(nextFree + lst.size());

    forAll(lst, elemI)
    {
        this->operator[](nextFree++) = lst[elemI];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::operator=(const T& t)
{
    UList<T>::operator=(t);
}


// ************************************************************************* //
