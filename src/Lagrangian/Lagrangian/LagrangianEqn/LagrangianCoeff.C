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

#include "LagrangianCoeff.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::initialise(const dimensionSet& dims)
{
    if (valid()) return;

    S_.set
    (
        new LagrangianSubField<Type>
        (
            IOobject
            (
                eqn_.name() + ":S" + (Implicit ? 'p' : 'u'),
                eqn_.mesh().mesh().time().name(),
                eqn_.mesh().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            eqn_.mesh(),
            dimensioned<Type>(dims, Zero)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, bool Implicit>
Foam::LagrangianCoeff<Type, Implicit>::LagrangianCoeff
(
    const LagrangianEqnBase& eqn
)
:
    eqn_(eqn)
{}


template<class Type, bool Implicit>
Foam::LagrangianCoeff<Type, Implicit>::LagrangianCoeff
(
    const LagrangianCoeff<Type, Implicit>& coeff
)
:
    tmp<LagrangianCoeff<Type, Implicit>>::refCount(),
    eqn_(coeff.eqn_),
    S_(coeff.S_, false)
{}


template<class Type, bool Implicit>
Foam::LagrangianCoeff<Type, Implicit>::LagrangianCoeff
(
    LagrangianCoeff<Type, Implicit>& coeff,
    const bool reuse
)
:
    eqn_(coeff.eqn_),
    S_(coeff.S_, reuse)
{}


template<class Type, bool Implicit>
Foam::LagrangianCoeff<Type, Implicit>::LagrangianCoeff
(
    LagrangianCoeff<Type, Implicit>&& coeff
)
:
    tmp<LagrangianCoeff<Type, Implicit>>::refCount(),
    eqn_(coeff.eqn_),
    S_(coeff.S_, true)
{}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
Foam::LagrangianCoeff<Type, Implicit>::LagrangianCoeff
(
    const LagrangianEqnBase& eqn,
    const LagrangianSubField<Type, PrimitiveField>& S
)
:
    eqn_(eqn),
    S_(S.clone())
{}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
Foam::LagrangianCoeff<Type, Implicit>::LagrangianCoeff
(
    const LagrangianEqnBase& eqn,
    const tmp<LagrangianSubField<Type, PrimitiveField>>& tS
)
:
    eqn_(eqn),
    S_(tS.isTmp() ? tS.ptr() : tS().clone().ptr())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, bool Implicit>
const Foam::LagrangianEqnBase&
Foam::LagrangianCoeff<Type, Implicit>::eqn() const
{
    return eqn_;
}


template<class Type, bool Implicit>
bool Foam::LagrangianCoeff<Type, Implicit>::valid() const
{
    return S_.valid();
}


template<class Type, bool Implicit>
const Foam::LagrangianSubField<Type>&
Foam::LagrangianCoeff<Type, Implicit>::S() const
{
    return S_();
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::negate()
{
    if (valid()) S_() = -S_();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator+=
(
    const LagrangianSubField<Type, PrimitiveField>& S
)
{
    initialise(S.dimensions());
    S_() += S;
}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator+=
(
    const tmp<LagrangianSubField<Type, PrimitiveField>>& tS
)
{
    operator+=(tS());
    tS.clear();
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator+=
(
    const LagrangianCoeff& coeff
)
{
    if (!coeff.valid()) return;
    operator+=(coeff.S_());
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator+=
(
    const dimensioned<Type>& dt
)
{
    initialise(dt.dimensions());
    S_() += dt;
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator+=(const zero)
{}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator-=
(
    const LagrangianSubField<Type, PrimitiveField>& S
)
{
    initialise(S.dimensions());
    S_() -= S;
}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator-=
(
    const tmp<LagrangianSubField<Type, PrimitiveField>>& tS
)
{
    operator-=(tS());
    tS.clear();
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator-=
(
    const LagrangianCoeff& coeff
)
{
    if (!coeff.valid()) return;
    operator-=(coeff.S_());
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator-=
(
    const dimensioned<Type>& dt
)
{
    initialise(dt.dimensions());
    S_() -= dt;
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator-=(const zero)
{}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator*=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    if (!valid()) return;
    S_() *= S;
}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator*=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    if (!valid()) return;
    S_() *= tS();
    tS.clear();
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator*=
(
    const dimensioned<scalar>& dt
)
{
    if (!valid()) return;
    S_() *= dt;
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator*=(const zero&)
{
    S_.clear();
}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator/=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    if (!valid()) return;
    S_() /= S;
}


template<class Type, bool Implicit>
template<template<class> class PrimitiveField>
void Foam::LagrangianCoeff<Type, Implicit>::operator/=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    if (!valid()) return;
    S_() /= tS();
    tS.clear();
}


template<class Type, bool Implicit>
void Foam::LagrangianCoeff<Type, Implicit>::operator/=
(
    const dimensioned<scalar>& dt
)
{
    if (!valid()) return;
    S_() /= dt;
}


// ************************************************************************* //
