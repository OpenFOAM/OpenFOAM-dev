/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class TypeB>
void Foam::LagrangianEqn<Type>::preOpCheck
(
    const tmp<LagrangianEqn<Type>>& tA,
    const tmp<LagrangianEqn<TypeB>>& tB
) const
{
    // Check that the meshes are the same
    if (&tA().mesh() != &tB().mesh())
    {
        FatalErrorInFunction
            << "Operating equations with different meshes"
            << exit(FatalError);
    }

    // Check that the referred fields are the same
    const regIOobject& aPsi =
        notNull(tA().psiSubSubRef_)
      ? static_cast<const regIOobject&>(tA().psiSubSubRef_)
      : notNull(tA().psiSubSub_)
      ? static_cast<const regIOobject&>(tA().psiSubSub_)
      : notNull(tA().psiSub_)
      ? static_cast<const regIOobject&>(tA().psiSub_)
      : NullObjectRef<regIOobject>();
    const regIOobject& bPsi =
        notNull(tB().psiSubSubRef_)
      ? static_cast<const regIOobject&>(tB().psiSubSubRef_)
      : notNull(tB().psiSubSub_)
      ? static_cast<const regIOobject&>(tB().psiSubSub_)
      : notNull(tB().psiSub_)
      ? static_cast<const regIOobject&>(tB().psiSub_)
      : NullObjectRef<regIOobject>();
    if (notNull(aPsi) && notNull(bPsi) && &aPsi != &bPsi)
    {
        FatalErrorInFunction
            << "Operating equations with different fields"
            << exit(FatalError);
    }
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::opDeltaT
(
    const tmp<LagrangianEqn<OtherType>>& tOther
)
{
    if (tOther().tDeltaT_.valid())
    {
        if (!tDeltaT_.valid())
        {
            tDeltaT_ =
                tmp<LagrangianSubScalarField>
                (
                    tOther().tDeltaT_,
                    tOther.isTmp()
                );
        }
        else if (&tDeltaT_() != &tOther().tDeltaT_())
        {
            FatalErrorInFunction
                << "Operating equations with different time-step fields"
                << exit(FatalError);
        }
    }
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::opClear
(
    LagrangianEqn<OtherType>& other
)
{
    // Remove old-time pointers so that they are not operated on on destruction
    other.deltaTSp0Ptr_ = NullObjectPtr<LagrangianDynamicField<scalar>>();
    other.S0Ptr_ = NullObjectPtr<LagrangianDynamicField<OtherType>>();
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::opClear
(
    const tmp<LagrangianEqn<OtherType>>& tOther
)
{
    if (tOther.isTmp())
    {
        opClear(tOther.ref());
    }
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::opFinalise
(
    LagrangianEqn<OtherType>& other
)
{
    // Store previous implicit time coefficient and source fields
    if (notNull(other.deltaTSp0Ptr_))
    {
        LagrangianSubSubField<scalar> deltaTSp0
        (
            other.mesh().sub(other.deltaTSp0Ptr_->internalField())
        );

        deltaTSp0 = Zero;

        if (other.deltaTSp.valid()) deltaTSp0 = other.deltaTSp.S();
    }

    if (notNull(other.S0Ptr_))
    {
        LagrangianSubSubField<Type> S0
        (
            other.mesh().sub(other.S0Ptr_->internalField())
        );

        other.mesh().group() == LagrangianGroup::none ? S0 = Zero : S0 *= -1;

        if (other.Su.valid()) S0 += other.Su.S();
        if (other.Sp.valid()) S0 += other.Sp.Su(other)->S();
    }

    opClear(other);
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::opFinalise
(
    const tmp<LagrangianEqn<OtherType>>& tOther
)
{
    if (tOther.isTmp())
    {
        opFinalise(tOther.ref());
    }
}


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
    psiSubSubRef_(NullObjectNonConstRef<LagrangianSubSubField<Type>>()),
    psiSubSub_(PsiRef<SubField>::ref(psi)),
    psiSub_(PsiRef<Field>::ref(psi)),
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
    psiSubSubRef_(NullObjectNonConstRef<LagrangianSubSubField<Type>>()),
    psiSubSub_(PsiRef<SubField>::ref(psi)),
    psiSub_(PsiRef<Field>::ref(psi)),
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
    psiSubSubRef_(psi),
    psiSubSub_(psi),
    psiSub_(NullObjectRef<LagrangianSubField<Type>>()),
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
Foam::LagrangianEqn<Type>::LagrangianEqn(const LagrangianSubMesh& mesh)
:
    LagrangianEqnBase(word::null, mesh),
    tDeltaT_(tmp<LagrangianSubScalarField>()),
    psiSubSubRef_(NullObjectNonConstRef<LagrangianSubSubField<Type>>()),
    psiSubSub_(NullObjectRef<LagrangianSubSubField<Type>>()),
    psiSub_(NullObjectRef<LagrangianSubField<Type>>()),
    deltaTSp0Ptr_(&NullObjectNonConstRef<LagrangianDynamicField<scalar>>()),
    S0Ptr_(&NullObjectNonConstRef<LagrangianDynamicField<Type>>()),
    deltaTSu(*this),
    deltaTSp(*this),
    Su(*this),
    Sp(*this)
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn(const LagrangianEqn<Type>& eqn)
:
    LagrangianEqn(tmp<LagrangianEqn<Type>>(eqn))
{}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn(const tmp<LagrangianEqn<Type>>& tEqn)
:
    LagrangianEqnBase(tEqn().name(), tEqn().mesh()),
    tDeltaT_
    (
        tEqn().tDeltaT_.valid()
      ? tmp<LagrangianSubScalarField>(tEqn().tDeltaT_, tEqn.isTmp())
      : tmp<LagrangianSubScalarField>()
    ),
    psiSubSubRef_(tEqn().psiSubSubRef_),
    psiSubSub_(tEqn().psiSubSub_),
    psiSub_(tEqn().psiSub_),
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
      ? LagrangianCoeff<Type, false>(*this, tEqn.ref().deltaTSu, true)
      : LagrangianCoeff<Type, false>(*this, tEqn().deltaTSu)
    ),
    deltaTSp
    (
        tEqn.isTmp()
      ? LagrangianCoeff<scalar, true>(*this, tEqn.ref().deltaTSp, true)
      : LagrangianCoeff<scalar, true>(*this, tEqn().deltaTSp)
    ),
    Su
    (
        tEqn.isTmp()
      ? LagrangianCoeff<Type, false>(*this, tEqn.ref().Su, true)
      : LagrangianCoeff<Type, false>(*this, tEqn().Su)
    ),
    Sp
    (
        tEqn.isTmp()
      ? LagrangianSp<Type>(*this, tEqn.ref().Sp, true)
      : LagrangianSp<Type>(*this, tEqn().Sp)
    )
{
    opClear(tEqn);
    tEqn.clear();
}


template<class Type>
Foam::LagrangianEqn<Type>::LagrangianEqn(LagrangianEqn<Type>&& eqn)
:
    LagrangianEqnBase(eqn.name(), eqn.mesh()),
    tDeltaT_(move(eqn.tDeltaT_)),
    psiSubSubRef_(eqn.psiSubSubRef_),
    psiSubSub_(eqn.psiSubSub_),
    psiSub_(eqn.psiSub_),
    deltaTSp0Ptr_(eqn.deltaTSp0Ptr_),
    S0Ptr_(eqn.S0Ptr_),
    deltaTSu(move(eqn.deltaTSu)),
    deltaTSp(move(eqn.deltaTSp)),
    Su(move(eqn.Su)),
    Sp(move(eqn.Sp))
{
    opClear(eqn);
}


template<class Type>
template<class EqOp>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const tmp<LagrangianEqn<Type>>& tA,
    const tmp<LagrangianEqn<Type>>& tB,
    const EqOp& eqOp
)
:
    LagrangianEqnBase
    (
        notNull(tA().psiSubSubRef_) && isNull(tB().psiSubSubRef_) ? tA().name()
      : isNull(tA().psiSubSubRef_) && notNull(tB().psiSubSubRef_) ? tB().name()
      : tA().name() == word::null ? tB().name()
      : tB().name() == word::null ? tA().name()
      : tA().name() == tB().name() ? tA().name()
      : word::null,
        tA().mesh()
    ),
    tDeltaT_
    (
        tA().tDeltaT_.valid()
      ? tmp<LagrangianSubScalarField>(tA().tDeltaT_, tA.isTmp())
      : tmp<LagrangianSubScalarField>()
    ),
    psiSubSubRef_
    (
        notNull(tA().psiSubSubRef_) ? tA().psiSubSubRef_
      : notNull(tB().psiSubSubRef_) ? tB().psiSubSubRef_
      : NullObjectNonConstRef<LagrangianSubSubField<Type>>()
    ),
    psiSubSub_
    (
        notNull(tA().psiSubSub_) ? tA().psiSubSub_
      : notNull(tB().psiSubSub_) ? tB().psiSubSub_
      : NullObjectRef<LagrangianSubSubField<Type>>()
    ),
    psiSub_
    (
        notNull(tA().psiSub_) ? tA().psiSub_
      : notNull(tB().psiSub_) ? tB().psiSub_
      : NullObjectRef<LagrangianSubField<Type>>()
    ),
    deltaTSp0Ptr_
    (
        tA.isTmp()
      ? tA().deltaTSp0Ptr_
      : NullObjectPtr<LagrangianDynamicField<scalar>>()
    ),
    S0Ptr_
    (
        tA.isTmp()
      ? tA().S0Ptr_
      : NullObjectPtr<LagrangianDynamicField<Type>>()
    ),
    deltaTSu
    (
        tA.isTmp()
      ? LagrangianCoeff<Type, false>(*this, tA.ref().deltaTSu, true)
      : LagrangianCoeff<Type, false>(*this, tA().deltaTSu)
    ),
    deltaTSp
    (
        tA.isTmp()
      ? LagrangianCoeff<scalar, true>(*this, tA.ref().deltaTSp, true)
      : LagrangianCoeff<scalar, true>(*this, tA().deltaTSp)
    ),
    Su
    (
        tA.isTmp()
      ? LagrangianCoeff<Type, false>(*this, tA.ref().Su, true)
      : LagrangianCoeff<Type, false>(*this, tA().Su)
    ),
    Sp
    (
        tA.isTmp()
      ? LagrangianSp<Type>(*this, tA.ref().Sp, true)
      : LagrangianSp<Type>(*this, tA().Sp)
    )
{
    preOpCheck(tA, tB);

    opDeltaT(tB);

    opClear(tA);
    tA.clear();

    // Combine the coefficients
    eqOp(*this, tB());

    // Copy/check the old-time pointers from/against the second equation
    if (tB.isTmp())
    {
        if (notNull(tB().deltaTSp0Ptr_))
        {
            if (notNull(deltaTSp0Ptr_) && deltaTSp0Ptr_ != tB().deltaTSp0Ptr_)
            {
                FatalErrorInFunction
                    << "Operating equations with different previous "
                    << "implicit time coefficients" << exit(FatalError);
            }

            if (isNull(deltaTSp0Ptr_)) deltaTSp0Ptr_ = tB.ref().deltaTSp0Ptr_;
        }

        if (notNull(tB().S0Ptr_))
        {
            if (notNull(S0Ptr_) && S0Ptr_ != tB().S0Ptr_)
            {
                FatalErrorInFunction
                    << "Operating equations with different previous "
                    << "sources" << exit(FatalError);
            }

            if (isNull(S0Ptr_)) S0Ptr_ = tB.ref().S0Ptr_;
        }
    }

    opClear(tB);
    tB.clear();
}


template<class Type>
template<class EqOp, class TypeB>
Foam::LagrangianEqn<Type>::LagrangianEqn
(
    const tmp<LagrangianEqn<Type>>& tA,
    const tmp<LagrangianEqn<TypeB>>& tB,
    const EqOp& eqOp
)
:
    LagrangianEqn(tA)
{
    preOpCheck(tmp<LagrangianEqn<Type>>(*this), tB);

    opDeltaT(tB);

    // Ignore the old-time pointers in the second equation. If they are of a
    // different type then they cannot relate to this equation.

    // Combine the coefficients
    eqOp(*this, tB());

    opClear(tB);
    tB.clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianEqn<Type>::~LagrangianEqn()
{
    opFinalise(*this);
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
    if (notNull(psiSubSub_)) return psiSubSub_;

    if (notNull(psiSub_)) return toSubField(psiSub_);

    FatalErrorInFunction
        << "Requested field from field-less equation"
        << exit(FatalError);

    return tmp<LagrangianSubSubField<Type>>();
}


template<class Type>
const Foam::word& Foam::LagrangianEqn<Type>::psiName() const
{
    if (notNull(psiSubSub_)) return psiSubSub_.name();

    if (notNull(psiSub_)) return psiSub_.name();

    return word::null;
}


template<class Type>
Foam::word Foam::LagrangianEqn<Type>::psiGroup() const
{
    if (notNull(psiSubSub_)) return psiSubSub_.group();

    if (notNull(psiSub_)) return psiSub_.group();

    return word::null;
}


template<class Type>
template<template<class> class PrimitiveField>
bool Foam::LagrangianEqn<Type>::isPsi
(
    const LagrangianSubField<Type, PrimitiveField>& otherPsi
) const
{
    const regIOobject& otherPsiIo = otherPsi;

    if (notNull(psiSubSub_))
    {
        return &otherPsiIo == &static_cast<const regIOobject&>(psiSubSub_);
    }

    if (notNull(psiSub_))
    {
        return &otherPsiIo == &static_cast<const regIOobject&>(psiSub_);
    }

    return false;
}


template<class Type>
bool Foam::LagrangianEqn<Type>::valid() const
{
    return deltaTSu.valid() || deltaTSp.valid() || Su.valid() || Sp.valid();
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianEqn<Type>::allSu() const
{
    if (Su.valid() && !tDeltaT_.valid())
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
    if (tDeltaT_.valid()) tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSu;

    return tResult;
}


template<class Type>
Foam::tmp<Foam::LagrangianSp<Type>> Foam::LagrangianEqn<Type>::allSp() const
{
    if (Sp.valid() && !tDeltaT_.valid())
    {
        FatalErrorInFunction
            << "Cannot evaluate equation " << name_ << " for field "
            << psiName() << " without a time-step" << exit(FatalError);
    }

    tmp<LagrangianSp<Type>> tResult(new LagrangianSp<Type>(*this));

    tResult.ref() += Sp;
    if (tDeltaT_.valid()) tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSp;

    return tResult;
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianEqn<Type>::allDiagonalSu() const
{
    tmp<LagrangianCoeff<Type, false>> tSpH(Sp.H());

    if ((Su.valid() || tSpH().valid()) && !tDeltaT_.valid())
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
    tResult.ref() += tSpH;
    if (tDeltaT_.valid()) tResult.ref() *= tDeltaT_();
    tResult.ref() += deltaTSu;

    return tResult;
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Foam::scalar, true>>
Foam::LagrangianEqn<Type>::allDiagonalSp() const
{
    tmp<LagrangianCoeff<scalar, true>> tSpA(Sp.A());

    if (tSpA().valid() && !tDeltaT_.valid())
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
    if (tDeltaT_.valid()) tResult.ref() *= tDeltaT_();
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

    if (isNull(psiSubSubRef_))
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
        psiSubSubRef_ = Zero;
        return;
    }

    psiSubSubRef_ = - (allSu()()/allSp()());

    if (!final)
    {
        deltaTSp0Ptr_ = NullObjectPtr<LagrangianDynamicField<scalar>>();
        S0Ptr_ = NullObjectPtr<LagrangianDynamicField<Type>>();
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::operator+=
(
    const LagrangianEqn<OtherType>& other
)
{
    deltaTSu += other.deltaTSu;
    deltaTSp += other.deltaTSp;
    Su += other.Su;
    Sp += other.Sp;
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::operator+=
(
    const tmp<LagrangianEqn<OtherType>>& tOther
)
{
    operator+=(tOther());
    tOther.clear();
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::operator-=
(
    const LagrangianEqn<OtherType>& other
)
{
    deltaTSu -= other.deltaTSu;
    deltaTSp -= other.deltaTSp;
    Su -= other.Su;
    Sp -= other.Sp;
}


template<class Type>
template<class OtherType>
void Foam::LagrangianEqn<Type>::operator-=
(
    const tmp<LagrangianEqn<OtherType>>& tOther
)
{
    operator-=(tOther());
    tOther.clear();
}


template<class Type>
template<template<class> class PrimitiveField>
void Foam::LagrangianEqn<Type>::operator*=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    opFinalise(*this);
    deltaTSu *= S;
    deltaTSp *= S;
    Su *= S;
    Sp *= S;
}


template<class Type>
template<template<class> class PrimitiveField>
void Foam::LagrangianEqn<Type>::operator*=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    this->operator*=(tS());
    tS.clear();
}


template<class Type>
void Foam::LagrangianEqn<Type>::operator*=(const dimensioned<scalar>& dt)
{
    opFinalise(*this);
    deltaTSu *= dt;
    deltaTSp *= dt;
    Su *= dt;
    Sp *= dt;
}


template<class Type>
void Foam::LagrangianEqn<Type>::operator*=(const zero& z)
{
    opFinalise(*this);
    deltaTSu *= z;
    deltaTSp *= z;
    Su *= z;
    Sp *= z;
}


template<class Type>
template<template<class> class PrimitiveField>
void Foam::LagrangianEqn<Type>::operator/=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    opFinalise(*this);
    deltaTSu /= S;
    deltaTSp /= S;
    Su /= S;
    Sp /= S;
}


template<class Type>
template<template<class> class PrimitiveField>
void Foam::LagrangianEqn<Type>::operator/=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    this->operator/=(tS());
    tS.clear();
}


template<class Type>
void Foam::LagrangianEqn<Type>::operator/=(const dimensioned<scalar>& dt)
{
    opFinalise(*this);
    deltaTSu /= dt;
    deltaTSp /= dt;
    Su /= dt;
    Sp /= dt;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator-
(
    const LagrangianEqn<Type>& eqn
)
{
    return - tmp<LagrangianEqn<Type>>(eqn);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator-
(
    const tmp<LagrangianEqn<Type>>& tEqn
)
{
    tmp<LagrangianEqn<Type>> tResult(tEqn, true);
    tResult.ref().deltaTSu.negate();
    tResult.ref().deltaTSp.negate();
    tResult.ref().Su.negate();
    tResult.ref().Sp.negate();
    return tResult;
}


#define LAGRANGIAN_EQN_EQN_OPERATOR(Op, EqOp)                                  \
                                                                               \
    template<class Type, class TypeB>                                          \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const LagrangianEqn<TypeB>& b                                          \
    )                                                                          \
    {                                                                          \
        return tmp<LagrangianEqn<Type>>(a) Op tmp<LagrangianEqn<TypeB>>(b);    \
    };                                                                         \
                                                                               \
    template<class Type, class TypeB>                                          \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tA,                                    \
        const LagrangianEqn<TypeB>& b                                          \
    )                                                                          \
    {                                                                          \
        return tA Op tmp<LagrangianEqn<TypeB>>(b);                             \
    };                                                                         \
                                                                               \
    template<class Type, class TypeB>                                          \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tA,                                    \
        const tmp<LagrangianEqn<TypeB>>& tB                                    \
    )                                                                          \
    {                                                                          \
        return                                                                 \
            tmp<LagrangianEqn<Type>>                                           \
            (                                                                  \
                new LagrangianEqn<Type>                                        \
                (                                                              \
                    tA,                                                        \
                    tB,                                                        \
                    []                                                         \
                    (                                                          \
                        LagrangianEqn<Type>& a,                                \
                        const LagrangianEqn<TypeB>& b                          \
                    )                                                          \
                    {                                                          \
                        a EqOp b;                                              \
                    }                                                          \
                )                                                              \
            );                                                                 \
    }                                                                          \
                                                                               \
    template<class Type, class TypeB>                                          \
    Foam::tmp<Foam::LagrangianEqn<Type>> Foam::operator Op                     \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const tmp<LagrangianEqn<TypeB>>& tB                                    \
    )                                                                          \
    {                                                                          \
        return tmp<LagrangianEqn<Type>>(a) Op tB;                              \
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

LAGRANGIAN_COMMUTATIVE_EQN_EQN_OPERATOR(+, +=)

LAGRANGIAN_EQN_EQN_OPERATOR(-, -=)

LAGRANGIAN_EQN_EQN_OPERATOR(==, -=)

#undef LAGRANGIAN_EQN_EQN_OPERATOR
#undef LAGRANGIAN_COMMUTATIVE_EQN_EQN_OPERATOR


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
        tmp<LagrangianEqn<Type>> tResult(tEqn, true);                          \
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
        tmp<LagrangianEqn<Type>> tResult(tEqn, true);                          \
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
        tResult.ref() EqOp field;                                              \
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
        tmp<LagrangianEqn<Type>> tResult(tEqn, true);                          \
        tResult.ref() EqOp field;                                              \
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
        tmp<LagrangianEqn<Type>> tResult(tEqn, true);                          \
        tResult.ref() EqOp tField;                                             \
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
        tResult.ref() EqOp tField;                                             \
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
