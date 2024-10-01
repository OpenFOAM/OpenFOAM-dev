/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "SubField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubList<Type>& list
)
:
    SubList<Type>(list)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list
)
:
    SubList<Type>(list, list.size())
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label subSize
)
:
    SubList<Type>(list, subSize)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label subSize,
    const label startIndex
)
:
    SubList<Type>(list, subSize, startIndex)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubField<Type>& sfield
)
:
    tmp<SubField<Type>>::refCount(),
    SubList<Type>(sfield)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    SubField<Type>& sfield,
    bool reuse
)
:
    tmp<SubField<Type>>::refCount(),
    SubList<Type>(sfield)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::Field<typename Foam::SubField<Type>::cmptType>>
Foam::SubField<Type>::component
(
    const direction d
) const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).component(d);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::SubField<Type>::T() const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).T();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SubField<Type>::operator=(const SubField<Type>& rhs)
{
    SubList<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const UList<Type>& rhs)
{
    SubList<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const tmp<Field<Type>>& rhs)
{
    SubList<Type>::operator=(rhs());
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const Type& rhs)
{
    SubList<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const zero)
{
    SubList<Type>::operator=(Zero);
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type>                                                           \
void Foam::SubField<Type>::operator op(const UList<TYPE>& f)                   \
{                                                                              \
    TFOR_ALL_F_OP_F(Type, *this, op, TYPE, f)                                  \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::SubField<Type>::operator op(const tmp<Field<TYPE>>& tf)             \
{                                                                              \
    operator op(tf());                                                         \
    tf.clear();                                                                \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::SubField<Type>::operator op(const TYPE& t)                          \
{                                                                              \
    TFOR_ALL_F_OP_S(Type, *this, op, TYPE, t)                                  \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


template<class Type>
template<class Form, Foam::direction Ncmpts>
inline void Foam::SubField<Type>::operator=
(
    const VectorSpace<Form, Type, Ncmpts>& rhs
)
{
    forAll(rhs, i)
    {
        this->operator[](i) = rhs[i];
    }
}


// ************************************************************************* //
