/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "oneOrTmp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::oneOrTmp<Type>::oneOrTmp()
:
    isOne_(true),
    tt_(nullptr)
{}


template<class Type>
inline Foam::oneOrTmp<Type>::oneOrTmp(const Type& t)
:
    oneOrTmp(tmp<Type>(t))
{}


template<class Type>
inline Foam::oneOrTmp<Type>::oneOrTmp(const tmp<Type>& tt)
:
    isOne_(false),
    tt_(tt)
{
    tt.clear();
}


template<class Type>
inline Foam::oneOrTmp<Type>::oneOrTmp(const oneOrTmp<Type>& oot)
:
    isOne_(oot.isOne_),
    tt_(oot.tt_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::oneOrTmp<Type>::~oneOrTmp()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
template<class Other>
void Foam::oneOrTmp<Type>::scale(Other& o) const
{
    if (!isOne_)
    {
        o *= tt_;
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::oneOrTmp<Type>::operator=(const one&)
{
    isOne_ = true;
    tt_.clear();
}


template<class Type>
void Foam::oneOrTmp<Type>::operator=(const tmp<Type>& tt)
{
    isOne_ = false;
    tt_ = tt;
}


template<class Type>
Foam::oneOrTmp<Type> Foam::oneOrTmp<Type>::operator()() const
{
    if (!isOne_)
    {
        return oneOrTmp<Type>(tt_());
    }
    else
    {
        return oneOrTmp();
    }
}


template<class Type>
template<class Other>
Foam::tmp<Other> Foam::oneOrTmp<Type>::operator*(const Other& o)
{
    if (!isOne_)
    {
        return tt_*o;
    }
    else
    {
        return tmp<Other>(o);
    }
}


template<class Type>
template<class Other>
Foam::tmp<Other> Foam::oneOrTmp<Type>::operator*(const tmp<Other>& to)
{
    if (!isOne_)
    {
        return tt_*to;
    }
    else
    {
        tmp<Other> tresult(to);

        to.clear();

        return tresult;
    }
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

template<class Type, class Other>
Foam::tmp<Other> Foam::operator*(const Other& a, const oneOrTmp<Type>& b)
{
    // Multiplication is commutative
    return b*a;
}


template<class Type, class Other>
Foam::tmp<Other> Foam::operator*(const tmp<Other>& ta, const oneOrTmp<Type>& b)
{
    // Multiplication is commutative
    return b*ta;
}


template<class Type, class Other>
void Foam::operator*=(Other& a, const oneOrTmp<Type>& b)
{
    b.scale(a);
}


// ************************************************************************* //
