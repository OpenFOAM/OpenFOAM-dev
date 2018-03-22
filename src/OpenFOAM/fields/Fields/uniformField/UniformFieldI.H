/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "UniformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::UniformField<Type>::UniformField(const Type& value)
:
    value_(value)
{}


template<class Type>
inline Foam::UniformField<Type>::operator Type() const
{
    return value_;
}


template<class Type>
inline Type Foam::UniformField<Type>::operator[](const label) const
{
    return value_;
}


template<class Type>
inline Foam::UniformField<Type> Foam::UniformField<Type>::field() const
{
    return UniformField(value_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline UniformField<Type> min
(
    const UniformField<Type>& u1,
    const UniformField<Type>& u2
)
{
    return UniformField<Type>(min(u1.operator Type(), u2.operator Type()));
}


template<class Type, class OtherType>
inline OtherType min(const UniformField<Type>& u, const OtherType& o)
{
    return min(u.operator Type(), o);
}


template<class Type, class OtherType>
inline OtherType min(const OtherType& o, const UniformField<Type>& u)
{
    return min(o, u.operator Type());
}


template<class Type>
inline UniformField<Type> max
(
    const UniformField<Type>& u1,
    const UniformField<Type>& u2
)
{
    return UniformField<Type>(max(u1.operator Type(), u2.operator Type()));
}


template<class Type, class OtherType>
inline OtherType max(const UniformField<Type>& u, const OtherType& o)
{
    return max(u.operator Type(), o);
}


template<class Type, class OtherType>
inline OtherType max(const OtherType& o, const UniformField<Type>& u)
{
    return max(o, u.operator Type());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
