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

\*---------------------------------------------------------------------------*/

#include "PtrField.H"
#include "expressionAssert.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PtrField<Type>::PtrField()
:
    PtrList<Type>()
{}


template<class Type>
Foam::PtrField<Type>::PtrField(const label size)
:
    PtrList<Type>(size)
{}


template<class Type>
Foam::PtrField<Type>::PtrField(const PtrField<Type>& f)
:
    tmp<PtrField<Type>>::refCount(),
    PtrList<Type>(f)
{}


template<class Type>
Foam::PtrField<Type>::PtrField(PtrField<Type>&& f)
:
    tmp<PtrField<Type>>::refCount(),
    PtrList<Type>(move(f))
{}


template<class Type>
Foam::PtrField<Type>::PtrField(const tmp<PtrField<Type>>& tf)
:
    PtrList<Type>(const_cast<PtrField<Type>&>(tf()), tf.isTmp())
{
    tf.clear();
}


template<class Type>
Foam::PtrField<Type>::PtrField(const PtrList<Type>& f)
:
    PtrList<Type>(f.size())
{
    forAll(*this, i)
    {
        // ***HGW Temporary clone of fv?PatchField for testing
        this->set(i, f[i].clone(f[i].internalField()));
    }
}


template<class Type>
template< class Expression, class ... Args, class>
Foam::PtrField<Type>::PtrField(const Expression& e, const Args& ... args)
:
    PtrList<Type>(expression::getFirst<expression::Size>(e))
{
    expression::assertSameAllContainerProperty<expression::Size>(e);

    forAll(*this, i)
    {
        PtrList<Type>::set
        (
            i,
            expression::New<Type, true>
            (
                expression::access(e, i),
                args ...
            ).ptr()
        );
    }
}


template<class Type>
Foam::tmp<Foam::PtrField<Type>> Foam::PtrField<Type>::clone() const
{
    return tmp<PtrField<Type>>(new PtrField<Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PtrField<Type>::negate()
{
    forAll(*this, i)
    {
        this->operator[](i).negate();
    }
}


template<class Type>
void Foam::PtrField<Type>::replace
(
    const direction d,
    const PtrField<cmptField>& sf
)
{
    forAll(*this, i)
    {
        this->operator[](i).replace(d, sf[i]);
    }
}


template<class Type>
void Foam::PtrField<Type>::replace
(
    const direction d,
    const cmptType& s
)
{
    forAll(*this, i)
    {
        this->operator[](i).replace(d, s);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::PtrField<Type>::operator=(const PtrField<Type>& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    forAll(*this, i)
    {
        this->operator[](i) = f[i];
    }
}


template<class Type>
void Foam::PtrField<Type>::operator=(PtrField<Type>&& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    PtrList<Type>::operator=(move(f));
}


template<class Type>
void Foam::PtrField<Type>::operator=(const tmp<PtrField>& tf)
{
    if (this == &(tf()))
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    PtrField* fieldPtr = tf.ptr();
    PtrList<Type>::transfer(*fieldPtr);
    delete fieldPtr;
}


template<class Type>
template<class Expression, class>
void Foam::PtrField<Type>::operator=(const Expression& e)
{
    expression::assertSameAllContainerProperty<expression::Size>
    (
        *this,
        e
    );

    forAll(*this, i)
    {
        this->operator[](i) = expression::access(e, i);
    }
}


template<class Type>
void Foam::PtrField<Type>::operator=
(
    const typename PtrField<Type>::pType& t
)
{
    forAll(*this, i)
    {
        this->operator[](i) = t;
    }
}


#define COMPUTED_ASSIGNMENT(TYPE, PTYPE, op)                                   \
                                                                               \
template<class Type>                                                           \
void Foam::PtrField<Type>::operator op(const PtrField<TYPE>& f)                \
{                                                                              \
    forAll(*this, i)                                                           \
    {                                                                          \
        this->operator[](i) op f[i];                                           \
    }                                                                          \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::PtrField<Type>::operator op                                         \
(                                                                              \
    const tmp<PtrField<TYPE>>& tf                                              \
)                                                                              \
{                                                                              \
    operator op(tf());                                                         \
    tf.clear();                                                                \
}                                                                              \
                                                                               \
template<class Type>                                                           \
template<class Expression, class>                                              \
void Foam::PtrField<Type>::operator op(const Expression& e)                    \
{                                                                              \
    /* Error if the expression is a different size */                          \
    expression::assertSameAllContainerProperty<expression::Size>               \
    (                                                                          \
        *this,                                                                 \
        e                                                                      \
    );                                                                         \
                                                                               \
    forAll(*this, i)                                                           \
    {                                                                          \
        this->operator[](i) op expression::access(e, i);                       \
    }                                                                          \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::PtrField<Type>::operator op(const PTYPE& t)                         \
{                                                                              \
    forAll(*this, i)                                                           \
    {                                                                          \
        this->operator[](i) op t;                                              \
    }                                                                          \
}


#define pType_ typename PtrField<Type>::pType
COMPUTED_ASSIGNMENT(Type, pType_, +=)
COMPUTED_ASSIGNMENT(Type, pType_, -=)
#undef pType_

#define cmptField_ typename PtrField<Type>::cmptField
#define pCmptType_ typename PtrField<Type>::pCmptType
COMPUTED_ASSIGNMENT(cmptField_, pCmptType_, *=)
COMPUTED_ASSIGNMENT(cmptField_, pCmptType_, /=)
#undef pCmptType_
#undef cmptField_

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const PtrField<Type>& f)
{
    os << static_cast<const PtrList<Type>&>(f);
    return os;
}


// ************************************************************************* //
