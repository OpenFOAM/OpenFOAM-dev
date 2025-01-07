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

#include "LagrangianEqn.H"
#include "GeometricField.H"
#include "toSubField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
template<template<class> class PrimitiveField>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const word& name,
    const tmp<LagrangianSubScalarField>& tDeltaT,
    const LagrangianSubField<Type, PrimitiveField>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqnBase(name, psi.mesh()),
    tDeltaT_(tDeltaT.valid() ? tDeltaT : tmp<LagrangianSubScalarField>()),
    psiSubSub_(PsiRef<SubField>::ref(psi)),
    psiSub_(PsiRef<Field>::ref(psi)),
    psiPtr_(NullObjectPtr<LagrangianSubSubField<Type>>()),
    deltaTSp0Ptr_(&deltaTSp0),
    S0Ptr_(&S0),
    deltaTSu(*this),
    deltaTSp(*this),
    Su(*this),
    Sp(*this)
{}


template<class Type>
template<template<class> class PrimitiveField>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const word& name,
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubField<Type, PrimitiveField>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqn
    (
        name,
        tmp<LagrangianSubScalarField>(deltaT),
        psi,
        deltaTSp0,
        S0
    )
{}


template<class Type>
template<template<class> class PrimitiveField>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const word& name,
    const LagrangianSubField<Type, PrimitiveField>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqn(name, tmp<LagrangianSubScalarField>(), psi, deltaTSp0, S0)
{}


template<class Type>
template<template<class> class PrimitiveField>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const tmp<LagrangianSubScalarField>& tDeltaT,
    const LagrangianSubField<Type, PrimitiveField>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqnBase(word::null, psi.mesh()),
    tDeltaT_(tDeltaT.valid() ? tDeltaT : tmp<LagrangianSubScalarField>()),
    psiSubSub_(PsiRef<SubField>::ref(psi)),
    psiSub_(PsiRef<Field>::ref(psi)),
    psiPtr_(NullObjectPtr<LagrangianSubSubField<Type>>()),
    deltaTSp0Ptr_(&deltaTSp0),
    S0Ptr_(&S0),
    deltaTSu(*this),
    deltaTSp(*this),
    Su(*this),
    Sp(*this)
{}


template<class Type>
template<template<class> class PrimitiveField>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubField<Type, PrimitiveField>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqn(tmp<LagrangianSubScalarField>(deltaT), psi, deltaTSp0, S0)
{}


template<class Type>
template<template<class> class PrimitiveField>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const LagrangianSubField<Type, PrimitiveField>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqn(tmp<LagrangianSubScalarField>(), psi, deltaTSp0, S0)
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const tmp<LagrangianSubScalarField>& tDeltaT,
    LagrangianSubSubField<Type>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqnBase(psi.name() + "Eqn", psi.mesh()),
    tDeltaT_(tDeltaT.valid() ? tDeltaT : tmp<LagrangianSubScalarField>()),
    psiSubSub_(psi),
    psiSub_(NullObjectRef<LagrangianSubField<Type>>()),
    psiPtr_(&psi),
    deltaTSp0Ptr_(&deltaTSp0),
    S0Ptr_(&S0),
    deltaTSu(*this),
    deltaTSp(*this),
    Su(*this),
    Sp(*this)
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const LagrangianSubScalarField& deltaT,
    LagrangianSubSubField<Type>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqn(tmp<LagrangianSubScalarField>(deltaT), psi, deltaTSp0, S0)
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    LagrangianSubSubField<Type>& psi,
    LagrangianDynamicField<scalar>& deltaTSp0,
    LagrangianDynamicField<Type>& S0
)
:
    LagrangianEqn(tmp<LagrangianSubScalarField>(), psi, deltaTSp0, S0)
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn(const LagrangianEqn<Type>& eqn)
:
    tmp<LagrangianEqn<Type>>::refCount(),
    LagrangianEqnBase(eqn),
    tDeltaT_(eqn.tDeltaT_),
    psiSubSub_(eqn.psiSubSub_),
    psiSub_(eqn.psiSub_),
    psiPtr_(eqn.psiPtr_),
    deltaTSp0Ptr_(NullObjectPtr<LagrangianDynamicField<scalar>>()),
    S0Ptr_(NullObjectPtr<LagrangianDynamicField<Type>>()),
    deltaTSu(eqn.deltaTSu),
    deltaTSp(eqn.deltaTSp),
    Su(eqn.Su),
    Sp(eqn.Sp)
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn(LagrangianEqn<Type>&& eqn)
:
    tmp<LagrangianEqn<Type>>::refCount(),
    LagrangianEqnBase(eqn),
    tDeltaT_(move(eqn.tDeltaT_)),
    psiSubSub_(eqn.psiSubSub_),
    psiSub_(eqn.psiSub_),
    psiPtr_(eqn.psiPtr_),
    deltaTSp0Ptr_(eqn.deltaTSp0Ptr_),
    S0Ptr_(eqn.S0Ptr_),
    deltaTSu(move(eqn.deltaTSu)),
    deltaTSp(move(eqn.deltaTSp)),
    Su(move(eqn.Su)),
    Sp(move(eqn.Sp))
{
    eqn.deltaTSp0Ptr_ = NullObjectPtr<LagrangianDynamicField<scalar>>();
    eqn.S0Ptr_ = NullObjectPtr<LagrangianDynamicField<Type>>();
}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn(const tmp<LagrangianEqn<Type>>& tEqn)
:
    LagrangianEqnBase(tEqn()),
    tDeltaT_(tEqn().tDeltaT_, tEqn.isTmp()),
    psiSubSub_(tEqn().psiSubSub_),
    psiSub_(tEqn().psiSub_),
    psiPtr_(tEqn().psiPtr_),
    deltaTSp0Ptr_
    (
        tEqn.isTmp()
      ? tEqn().deltaTSp0Ptr_
      : NullObjectPtr<LagrangianDynamicField<scalar>>()
    ),
    S0Ptr_
    (
        tEqn.isTmp()
      ? tEqn().S0Ptr_
      : NullObjectPtr<LagrangianDynamicField<Type>>()
    ),
    deltaTSu
    (
        tEqn.isTmp()
      ? LagrangianCoeff<Type, false>(tEqn.ref().deltaTSu, true)
      : LagrangianCoeff<Type, false>(tEqn().deltaTSu)
    ),
    deltaTSp
    (
        tEqn.isTmp()
      ? LagrangianCoeff<scalar, true>(tEqn.ref().deltaTSp, true)
      : LagrangianCoeff<scalar, true>(tEqn().deltaTSp)
    ),
    Su
    (
        tEqn.isTmp()
      ? LagrangianCoeff<Type, false>(tEqn.ref().Su, true)
      : LagrangianCoeff<Type, false>(tEqn().Su)
    ),
    Sp
    (
        tEqn.isTmp()
      ? LagrangianSp<Type>(tEqn.ref().Sp, true)
      : LagrangianSp<Type>(tEqn().Sp)
    )
{
    if (tEqn.isTmp())
    {
        tEqn.ref().deltaTSp0Ptr_ =
            NullObjectPtr<LagrangianDynamicField<scalar>>();
    }
    if (tEqn.isTmp())
    {
        tEqn.ref().S0Ptr_ =
            NullObjectPtr<LagrangianDynamicField<Type>>();
    }

    tEqn.clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianEqn<Type>::~LagrangianEqn()
{
    if (notNull(deltaTSp0Ptr_))
    {
        LagrangianSubSubField<scalar> deltaTSp0
        (
            mesh().sub(deltaTSp0Ptr_->internalField())
        );

        deltaTSp0 = Zero;

        if (deltaTSp.valid()) deltaTSp0 = deltaTSp.S();
    }

    if (notNull(S0Ptr_))
    {
        LagrangianSubSubField<Type> S0
        (
            mesh().sub(S0Ptr_->internalField())
        );

        mesh().group() == LagrangianGroup::none ? S0 = Zero : S0 *= -1;

        if (Su.valid()) S0 += Su.S();
        if (Sp.valid()) S0 += Sp.Su(*this)->S();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::LagrangianEqn<Type>::NewEmpty
(
    const LagrangianEqn<Type>& eqn
)
{
    return
        tmp<LagrangianEqn<Type>>
        (
            notNull(eqn.psiSubSub_)
          ? new LagrangianEqn<Type>(eqn.psiSubSub_)
          : new LagrangianEqn<Type>(eqn.psiSub_)
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::LagrangianEqn<Type>::psi() const
{
    return
        notNull(psiSubSub_)
      ? tmp<LagrangianSubSubField<Type>>(psiSubSub_)
      : toSubField(psiSub_);
}


template<class Type>
const Foam::regIOobject& Foam::LagrangianEqn<Type>::psiIo() const
{
    return
        notNull(psiSubSub_)
      ? static_cast<const regIOobject&>(psiSubSub_)
      : static_cast<const regIOobject&>(psiSub_);
}


template<class Type>
const Foam::word& Foam::LagrangianEqn<Type>::psiName() const
{
    return psiIo().name();
}


template<class Type>
template<template<class> class PrimitiveField>
bool Foam::LagrangianEqn<Type>::isPsi
(
    const LagrangianSubField<Type, PrimitiveField>& psi
) const
{
    return &static_cast<const regIOobject&>(psi) == &psiIo();
}


template<class Type>
void Foam::LagrangianEqn<Type>::op(const LagrangianEqn<Type>& other)
{
    if
    (
        tDeltaT_.valid()
     && other.tDeltaT_.valid()
     && &tDeltaT_() != &other.tDeltaT_()
    )
    {
        FatalErrorInFunction
            << "Combining equations with different time-step fields"
            << exit(FatalError);
    }

    if (&psiIo() != &other.psiIo())
    {
        FatalErrorInFunction
            << "Combining equations with different fields"
            << exit(FatalError);
    }

    if (name_ == word::null || (isNull(psiPtr_) && notNull(other.psiPtr_)))
    {
        name_ = other.name_;
    }
    else if (name_ != other.name_)
    {
        name_ = word::null;
    }

    if (!tDeltaT_.valid() && other.tDeltaT_.valid())
    {
        tDeltaT_ = tmp<LagrangianSubScalarField>(other.tDeltaT_);
    }

    if (isNull(psiPtr_) && notNull(other.psiPtr_))
    {
        psiPtr_ = other.psiPtr_;
    }
}


template<class Type>
void Foam::LagrangianEqn<Type>::setPrevious
(
    const tmp<LagrangianEqn<Type>>& tOther
)
{
    if (!tOther.isTmp()) return;

    if (notNull(tOther().deltaTSp0Ptr_))
    {
        if (notNull(deltaTSp0Ptr_) && deltaTSp0Ptr_ != tOther().deltaTSp0Ptr_)
        {
            FatalErrorInFunction
                << "Operating tmp equations with different previous "
                << "implicit time coefficients" << exit(FatalError);
        }

        if (isNull(deltaTSp0Ptr_)) deltaTSp0Ptr_ = tOther.ref().deltaTSp0Ptr_;

        tOther.ref().deltaTSp0Ptr_ =
            NullObjectPtr<LagrangianDynamicField<scalar>>();
    }

    if (notNull(tOther().S0Ptr_))
    {
        if (notNull(S0Ptr_) && S0Ptr_ != tOther().S0Ptr_)
        {
            FatalErrorInFunction
                << "Operating tmp equations with different previous "
                << "sources" << exit(FatalError);
        }

        if (isNull(S0Ptr_)) S0Ptr_ = tOther.ref().S0Ptr_;

        tOther.ref().S0Ptr_ =
            NullObjectPtr<LagrangianDynamicField<Type>>();
    }

    tOther.clear();
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianEqn<Type>::allSu() const
{
    if (!tDeltaT_.valid())
    {
        FatalErrorInFunction
            << "Cannot evaluate equation " << name_ << " for field "
            << psiName() << " without a time-step" << exit(FatalError);
    }

    tmp<LagrangianCoeff<Type, false>> tResult
    (
        new LagrangianCoeff<Type, false>(*this)
    );

    tResult.ref() += Su;
    tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSu;

    return tResult;
}


template<class Type>
Foam::tmp<Foam::LagrangianSp<Type>> Foam::LagrangianEqn<Type>::allSp() const
{
    if (!tDeltaT_.valid())
    {
        FatalErrorInFunction
            << "Cannot evaluate equation " << name_ << " for field "
            << psiName() << " without a time-step" << exit(FatalError);
    }

    tmp<LagrangianSp<Type>> tResult(new LagrangianSp<Type>(*this));

    tResult.ref() += Sp;
    tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSp;

    return tResult;
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianEqn<Type>::allDiagonalSu() const
{
    if (!tDeltaT_.valid())
    {
        FatalErrorInFunction
            << "Cannot evaluate equation " << name_ << " for field "
            << psiName() << " without a time-step" << exit(FatalError);
    }

    tmp<LagrangianCoeff<Type, false>> tResult
    (
        new LagrangianCoeff<Type, false>(*this)
    );

    tResult.ref() += Su;
    tResult.ref() += Sp.H();
    tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSu;

    return tResult;
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Foam::scalar, true>>
Foam::LagrangianEqn<Type>::allDiagonalSp() const
{
    if (!tDeltaT_.valid())
    {
        FatalErrorInFunction
            << "Cannot evaluate equation " << name_ << " for field "
            << psiName() << " without a time-step" << exit(FatalError);
    }

    tmp<LagrangianCoeff<scalar, true>> tResult
    (
        new LagrangianCoeff<scalar, true>(*this)
    );

    tResult.ref() += Sp.A();
    tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSp;

    return tResult;
}


template<class Type>
void Foam::LagrangianEqn<Type>::solve(const bool final)
{
    if (!tDeltaT_.valid())
    {
        FatalErrorInFunction
            << "Cannot solve equation " << name_ << " for field "
            << psiName() << " without a time-step" << exit(FatalError);
    }

    if (isNull(psiPtr_))
    {
        FatalErrorInFunction
            << "Cannot solve equation " << name_ << " for constant field "
            << psiName() << exit(FatalError);
    }

    if (!deltaTSp.valid() && !Sp.valid())
    {
        FatalErrorInFunction
            << "Cannot solve equation " << name_ << " for field "
            << psiName() << " with no diagonal" << exit(FatalError);
    }

    if (!deltaTSu.valid() && !Su.valid())
    {
        *psiPtr_ = Zero;
        return;
    }

    *psiPtr_ = - (allSu()()/allSp()());

    if (!final)
    {
        deltaTSp0Ptr_ = NullObjectPtr<LagrangianDynamicField<scalar>>();
        S0Ptr_ = NullObjectPtr<LagrangianDynamicField<Type>>();
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianEqn<Type>::operator+=(const LagrangianEqn<Type>& other)
{
    deltaTSu += other.deltaTSu;
    deltaTSp += other.deltaTSp;
    Su += other.Su;
    Sp += other.Sp;
}


template<class Type>
void Foam::LagrangianEqn<Type>::operator+=
(
    const tmp<LagrangianEqn<Type>>& tOther
)
{
    operator+=(tOther());
    tOther.clear();
}


template<class Type>
void Foam::LagrangianEqn<Type>::operator-=(const LagrangianEqn<Type>& other)
{
    deltaTSu -= other.deltaTSu;
    deltaTSp -= other.deltaTSp;
    Su -= other.Su;
    Sp -= other.Sp;
}


template<class Type>
void Foam::LagrangianEqn<Type>::operator-=
(
    const tmp<LagrangianEqn<Type>>& tOther
)
{
    operator-=(tOther());
    tOther.clear();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator-
(
    const LagrangianEqn<Type>& eqn
)
{
    tmp<LagrangianEqn<Type>> tResult(new LagrangianEqn<Type>(eqn));
    tResult.ref().deltaTSu.negate();
    tResult.ref().deltaTSp.negate();
    tResult.ref().Su.negate();
    tResult.ref().Sp.negate();
    return tResult.ref();
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator-
(
    const tmp<LagrangianEqn<Type>>& tEqn
)
{
    tmp<LagrangianEqn<Type>> tResult(-tEqn());
    tEqn.clear();
    return tResult;
}


#define LAGRANGIAN_EQN_EQN_OPERATOR(Op, EqOp)                                  \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const LagrangianEqn<Type>& b                                           \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(new LagrangianEqn<Type>(a));          \
        tResult.ref().op(b);                                                   \
        tResult.ref() EqOp b;                                                  \
        return tResult;                                                        \
    };                                                                         \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tA,                                    \
        const LagrangianEqn<Type>& b                                           \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(tA);                                  \
        tResult.ref().op(b);                                                   \
        tResult.ref() EqOp b;                                                  \
        return tResult;                                                        \
    };                                                                         \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tA,                                    \
        const tmp<LagrangianEqn<Type>>& tB                                     \
    )                                                                          \
    {                                                                          \
                                                                               \
        if (!tA.isTmp()) return tA() Op tB;                                    \
        if (!tB.isTmp()) return tA Op tB();                                    \
        tmp<LagrangianEqn<Type>> tResult(tA Op tB());                          \
        tResult.ref().setPrevious(tB);                                         \
        return tResult;                                                        \
    }

#define LAGRANGIAN_COMMUTATIVE_EQN_EQN_OPERATOR(Op, EqOp)                      \
                                                                               \
    LAGRANGIAN_EQN_EQN_OPERATOR(Op, EqOp)                                      \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const tmp<LagrangianEqn<Type>>& tB                                     \
    )                                                                          \
    {                                                                          \
        return tB Op a;                                                        \
    }

#define LAGRANGIAN_NON_COMMUTATIVE_EQN_EQN_OPERATOR(Op, EqOp)                  \
                                                                               \
    LAGRANGIAN_EQN_EQN_OPERATOR(Op, EqOp)                                      \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const tmp<LagrangianEqn<Type>>& tB                                     \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(a Op tB());                           \
        tResult.ref().setPrevious(tB);                                         \
        return tResult;                                                        \
    }

LAGRANGIAN_COMMUTATIVE_EQN_EQN_OPERATOR(+, +=)

LAGRANGIAN_NON_COMMUTATIVE_EQN_EQN_OPERATOR(-, -=)

LAGRANGIAN_NON_COMMUTATIVE_EQN_EQN_OPERATOR(==, -=)

#undef LAGRANGIAN_EQN_EQN_OPERATOR
#undef LAGRANGIAN_COMMUTATIVE_EQN_EQN_OPERATOR
#undef LAGRANGIAN_NON_COMMUTATIVE_EQN_EQN_OPERATOR


#define LAGRANGIAN_EQN_FIELD_OPERATOR(Op, EqOp, LagrangianSubField)            \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const LagrangianSubField<Type>& field                                  \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(new LagrangianEqn<Type>(eqn));        \
        tResult.ref().Su EqOp field;                                           \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const LagrangianSubField<Type>& field                                  \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(tEqn);                                \
        tResult.ref().Su EqOp field;                                           \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const tmp<LagrangianSubField<Type>>& tField                            \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(tEqn);                                \
        tResult.ref().Su EqOp tField;                                          \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const tmp<LagrangianSubField<Type>>& tField                            \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(new LagrangianEqn<Type>(eqn));        \
        tResult.ref().Su EqOp tField;                                          \
        return tResult;                                                        \
    }

#define LAGRANGIAN_FIELD_EQN_OPERATOR(Op, EqOp, LagrangianSubField)            \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianSubField<Type>& field,                                 \
        const LagrangianEqn<Type>& eqn                                         \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult =                                     \
            LagrangianEqn<Type>::NewEmpty(eqn);                                \
        tResult.ref() += field;                                                \
        tResult.ref() EqOp eqn;                                                \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianSubField<Type>& field,                                 \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult =                                     \
            LagrangianEqn<Type>::NewEmpty(tEqn());                             \
        tResult.ref() += field;                                                \
        tResult.ref() EqOp tEqn;                                               \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianSubField<Type>>& tField,                           \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult =                                     \
            LagrangianEqn<Type>::NewEmpty(tEqn());                             \
        tResult.ref() += tField;                                               \
        tResult.ref() EqOp tEqn;                                               \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianSubField<Type>>& tField,                           \
        const LagrangianEqn<Type>& eqn                                         \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult =                                     \
            LagrangianEqn<Type>::NewEmpty(eqn);                                \
        tResult.ref() += tField;                                               \
        tResult.ref() EqOp eqn;                                                \
        return tResult;                                                        \
    }

//- Addition
LAGRANGIAN_EQN_FIELD_OPERATOR(+, +=, LagrangianSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(+, +=, LagrangianSubField)
LAGRANGIAN_EQN_FIELD_OPERATOR(+, +=, LagrangianSubSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(+, +=, LagrangianSubSubField)

//- Subtraction
LAGRANGIAN_EQN_FIELD_OPERATOR(-, -=, LagrangianSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(-, -=, LagrangianSubField)
LAGRANGIAN_EQN_FIELD_OPERATOR(-, -=, LagrangianSubSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(-, -=, LagrangianSubSubField)

//- Set-equal-to
LAGRANGIAN_EQN_FIELD_OPERATOR(==, -=, LagrangianSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(==, -=, LagrangianSubField)
LAGRANGIAN_EQN_FIELD_OPERATOR(==, -=, LagrangianSubSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(==, -=, LagrangianSubSubField)

#undef LAGRANGIAN_EQN_FIELD_OPERATOR
#undef LAGRANGIAN_FIELD_EQN_OPERATOR


#define LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(Op, EqOp, LagrangianSubField)     \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const LagrangianSubField<scalar>& field                                \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(new LagrangianEqn<Type>(eqn));        \
        tResult.ref().deltaTSu EqOp field;                                     \
        tResult.ref().deltaTSp EqOp field;                                     \
        tResult.ref().Su EqOp field;                                           \
        tResult.ref().Sp EqOp field;                                           \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const LagrangianSubField<scalar>& field                                \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(tEqn);                                \
        tResult.ref().deltaTSu EqOp field;                                     \
        tResult.ref().deltaTSp EqOp field;                                     \
        tResult.ref().Su EqOp field;                                           \
        tResult.ref().Sp EqOp field;                                           \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const tmp<LagrangianSubField<scalar>>& tField                          \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(tEqn);                                \
        tResult.ref().deltaTSu EqOp tField();                                  \
        tResult.ref().deltaTSp EqOp tField();                                  \
        tResult.ref().Su EqOp tField();                                        \
        tResult.ref().Sp EqOp tField();                                        \
        tField.clear();                                                        \
        return tResult;                                                        \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const tmp<LagrangianSubField<scalar>>& tField                          \
    )                                                                          \
    {                                                                          \
        tmp<LagrangianEqn<Type>> tResult(new LagrangianEqn<Type>(eqn));        \
        tResult.ref().deltaTSu EqOp tField();                                  \
        tResult.ref().deltaTSp EqOp tField();                                  \
        tResult.ref().Su EqOp tField();                                        \
        tResult.ref().Sp EqOp tField();                                        \
        tField.clear();                                                        \
        return tResult;                                                        \
    }

#define LAGRANGIAN_COMMUTATIVE_SCALAR_FIELD_EQN_OPERATOR(Op,LagrangianSubField)\
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianSubField<scalar>& field,                               \
        const LagrangianEqn<Type>& eqn                                         \
    )                                                                          \
    {                                                                          \
        return eqn Op field;                                                   \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianSubField<scalar>& field,                               \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    )                                                                          \
    {                                                                          \
        return tEqn Op field;                                                  \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianSubField<scalar>>& tField,                         \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    )                                                                          \
    {                                                                          \
        return tEqn Op tField;                                                 \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianSubField<scalar>>& tField,                         \
        const LagrangianEqn<Type>& eqn                                         \
    )                                                                          \
    {                                                                          \
        return eqn Op tField;                                                  \
    }

//- Multiplication
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(*, *=, LagrangianSubField)
LAGRANGIAN_COMMUTATIVE_SCALAR_FIELD_EQN_OPERATOR(*, LagrangianSubField)
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(*, *=, LagrangianSubSubField)
LAGRANGIAN_COMMUTATIVE_SCALAR_FIELD_EQN_OPERATOR(*, LagrangianSubSubField)

//- Division
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(/, /=, LagrangianSubField)
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(/, /=, LagrangianSubSubField)

#undef LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR
#undef LAGRANGIAN_COMMUTATIVE_SCALAR_FIELD_EQN_OPERATOR

// ************************************************************************* //
