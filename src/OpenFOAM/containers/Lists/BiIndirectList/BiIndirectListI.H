/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
inline Foam::BiIndirectList<T>::BiIndirectList
(
    const UList<T>& posList,
    const UList<T>& negList,
    const labelUList& addr
)
:
    posList_(const_cast<UList<T>&>(posList)),
    negList_(const_cast<UList<T>&>(negList)),
    addressing_(addr)
{}


template<class T>
inline Foam::BiIndirectList<T>::BiIndirectList
(
    const UList<T>& posList,
    const UList<T>& negList,
    List<label>&& addr
)
:
    posList_(const_cast<UList<T>&>(posList)),
    negList_(const_cast<UList<T>&>(negList)),
    addressing_(move(addr))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::BiIndirectList<T>::size() const
{
    return addressing_.size();
}


template<class T>
inline bool Foam::BiIndirectList<T>::empty() const
{
    return addressing_.empty();
}


template<class T>
inline const Foam::UList<T>& Foam::BiIndirectList<T>::posList() const
{
    return posList_;
}


template<class T>
inline const Foam::UList<T>& Foam::BiIndirectList<T>::negList() const
{
    return negList_;
}


template<class T>
inline const Foam::List<Foam::label>& Foam::BiIndirectList<T>::addressing()
 const
{
    return addressing_;
}


template<class T>
inline void Foam::BiIndirectList<T>::resetAddressing
(
    const labelUList& addr
)
{
    addressing_ = addr;
}


template<class T>
inline Foam::label Foam::BiIndirectList<T>::posIndex(const label i)
{
    return i;
}


template<class T>
inline Foam::label Foam::BiIndirectList<T>::negIndex(const label i)
{
    return -i-1;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T> Foam::BiIndirectList<T>::operator()() const
{
    List<T> result(size());

    forAll(*this, i)
    {
        result[i] = operator[](i);
    }

    return result;
}


template<class T>
inline T& Foam::BiIndirectList<T>::operator[](const label i)
{
    label index = addressing_[i];

    if (index >= 0)
    {
        return posList_[index];
    }
    else
    {
        return negList_[-index-1];
    }
}


template<class T>
inline const T& Foam::BiIndirectList<T>::operator[](const label i) const
{
    label index = addressing_[i];

    if (index >= 0)
    {
        return posList_[index];
    }
    else
    {
        return negList_[-index-1];
    }
}


template<class T>
inline void Foam::BiIndirectList<T>::operator=(List<T>&& ae)
{
    if (addressing_.size() != ae.size())
    {
        FatalErrorInFunction
            << "Addressing and list of addressed elements "
               "have different sizes: "
            << addressing_.size() << " " << ae.size()
            << abort(FatalError);
    }

    addressing_ = move(ae);
}


template<class T>
inline void Foam::BiIndirectList<T>::operator=(const UList<T>& ae)
{
    if (addressing_.size() != ae.size())
    {
        FatalErrorInFunction
            << "Addressing and list of addressed elements "
               "have different sizes: "
            << addressing_.size() << " " << ae.size()
            << abort(FatalError);
    }

    forAll(addressing_, i)
    {
        operator[](i) = ae[i];
    }
}


template<class T>
inline void Foam::BiIndirectList<T>::operator=(const T& t)
{
    forAll(addressing_, i)
    {
        operator[](i) = t;
    }
}


// ************************************************************************* //
