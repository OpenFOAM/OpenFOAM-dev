/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "RemoteData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::RemoteData<Type>::RemoteData()
:
    remote(),
    data()
{}


template<class Type>
inline Foam::RemoteData<Type>::RemoteData
(
    const label p,
    const label e,
    const Type& d
)
:
    remote(p, e),
    data(d)
{}


template<class Type>
inline Foam::RemoteData<Type>::RemoteData(const label e, const Type& d)
:
    remote(e),
    data(d)
{}


template<class Type>
inline Foam::RemoteData<Type>::RemoteData(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline Foam::RemoteData<Type> Foam::RemoteData<Type>::firstProcOp::operator()
(
    const RemoteData<Type>& a,
    const RemoteData<Type>& b
) const
{
    if (b.proci != -1 && (a.proci == -1 || a.proci < b.proci))
    {
        return b;
    }
    else
    {
        return a;
    }
}


template<class Type>
inline void Foam::RemoteData<Type>::firstProcEqOp::operator()
(
    RemoteData<Type>& a,
    const RemoteData<Type>& b
) const
{
    if (b.proci != -1 && (a.proci == -1 || a.proci < b.proci))
    {
        a = b;
    }
}


template<class Type>
inline void Foam::RemoteData<Type>::smallestEqOp::operator()
(
    RemoteData<Type>& a,
    const RemoteData<Type>& b
) const
{
    if (b.proci != -1 && (a.proci == -1 || b.data < a.data))
    {
        a = b;
    }
}


template<class Type>
inline void Foam::RemoteData<Type>::greatestEqOp::operator()
(
    RemoteData<Type>& a,
    const RemoteData<Type>& b
) const
{
    if (b.proci != -1 && (a.proci == -1 || b.data > a.data))
    {
        a = b;
    }
}


template<class Type>
inline void Foam::RemoteData<Type>::smallestFirstEqOp::operator()
(
    RemoteData<Type>& a,
    const RemoteData<Type>& b
) const
{
    if (b.proci != -1 && (a.proci == -1 || b.data.first() < a.data.first()))
    {
        a = b;
    }
}


template<class Type>
inline void Foam::RemoteData<Type>::greatestFirstEqOp::operator()
(
    RemoteData<Type>& a,
    const RemoteData<Type>& b
) const
{
    if (b.proci != -1 && (a.proci == -1 || b.data.first() > a.data.first()))
    {
        a = b;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::operator==
(
    const RemoteData<Type>& a,
    const RemoteData<Type>& b
)
{
    return
        static_cast<const remote&>(a) == static_cast<const remote&>(b)
     && a.data == b.data;
}


template<class Type>
inline bool Foam::operator!=
(
    const RemoteData<Type>& a,
    const RemoteData<Type>& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Type>
inline Foam::Ostream& Foam::operator<<(Ostream& os, const RemoteData<Type>& p)
{
    return os << static_cast<const remote&>(p) << token::SPACE << p.data;
}


template<class Type>
inline Foam::Istream& Foam::operator>>(Istream& is, RemoteData<Type>& p)
{
    return is >> static_cast<remote&>(p) >> p.data;
}


// ************************************************************************* //
