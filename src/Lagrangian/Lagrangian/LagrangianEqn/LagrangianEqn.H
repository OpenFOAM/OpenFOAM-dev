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

Class
    Foam::LagrangianEqn

Description
    This class stores the coefficients of a Lagrangian equation, and
    facilitates solving that equation and updating the associated field. It is
    designed to behave and be used similarly to fvMatrix.

SourceFiles
    LagrangianEqn.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianEqn_H
#define LagrangianEqn_H

#include "LagrangianSp.H"
#include "LagrangianEqnBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LagrangianEqn Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LagrangianEqn
:
    public tmp<LagrangianEqn<Type>>::refCount,
    public LagrangianEqnBase
{
    // Private Data

        //- The time-step field
        tmp<LagrangianSubScalarField> tDeltaT_;

        //- Reference to the sub-sub-field. One of this and psiSub will be
        //  valid and the other will be null.
        const LagrangianSubSubField<Type>& psiSubSub_;

        //- Reference to the sub-field. One of this and psiSubSub will be
        //  valid and the other will be null.
        const LagrangianSubField<Type>& psiSub_;

        //- Non-const pointer to the field. Might be null.
        LagrangianSubSubField<Type>* psiPtr_;

        //- Previous implicit time coefficient field into which to cache the
        //  implicit time coefficient on destruction. Might be null.
        LagrangianDynamicField<scalar>* deltaTSp0Ptr_;

        //- Previous source field into which to cache the source on
        //  destruction. Might be null.
        LagrangianDynamicField<Type>* S0Ptr_;


    // Private Classes

        //- Helper class to create the psiSub/psiSubSub references
        template<template<class> class PrimitiveField>
        struct PsiRef;


public:

    // Public Data

        //- Explicit time-coefficient
        LagrangianCoeff<Type, false> deltaTSu;

        //- Implicit time-coefficient
        LagrangianCoeff<scalar, true> deltaTSp;

        //- Explicit coefficient
        LagrangianCoeff<Type, false> Su;

        //- Implicit coefficient
        LagrangianSp<Type> Sp;


    // Constructors

        //- Construct for a const field and a tmp time-step with a name
        template<template<class> class PrimitiveField>
        LagrangianEqn
        (
            const word& name,
            const tmp<LagrangianSubScalarField>& tDeltaT,
            const LagrangianSubField<Type, PrimitiveField>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a const field and a time-step with a name
        template<template<class> class PrimitiveField>
        LagrangianEqn
        (
            const word& name,
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubField<Type, PrimitiveField>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a const field with a name
        template<template<class> class PrimitiveField>
        LagrangianEqn
        (
            const word& name,
            const LagrangianSubField<Type, PrimitiveField>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a const field and a tmp time-step. Name will be null.
        template<template<class> class PrimitiveField>
        LagrangianEqn
        (
            const tmp<LagrangianSubScalarField>& tDeltaT,
            const LagrangianSubField<Type, PrimitiveField>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a const field and a time-step. Name will be null.
        template<template<class> class PrimitiveField>
        LagrangianEqn
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubField<Type, PrimitiveField>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a const field. Name will be null.
        template<template<class> class PrimitiveField>
        LagrangianEqn
        (
            const LagrangianSubField<Type, PrimitiveField>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a non-const field and a tmp time-step. Will be named
        //  automatically as the field name plus the suffix "Eqn".
        LagrangianEqn
        (
            const tmp<LagrangianSubScalarField>& tDeltaT,
            LagrangianSubSubField<Type>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a non-const field and a time-step. Will be named
        //  automatically as the field name plus the suffix "Eqn".
        LagrangianEqn
        (
            const LagrangianSubScalarField& deltaT,
            LagrangianSubSubField<Type>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Construct for a non-const field. Will be named
        //  automatically as the field name plus the suffix "Eqn".
        LagrangianEqn
        (
            LagrangianSubSubField<Type>& psi,
            LagrangianDynamicField<scalar>& deltaTSp0 =
                NullObjectNonConstRef<LagrangianDynamicField<scalar>>(),
            LagrangianDynamicField<Type>& S0 =
                NullObjectNonConstRef<LagrangianDynamicField<Type>>()
        );

        //- Copy construct
        LagrangianEqn(const LagrangianEqn<Type>& eqn);

        //- Move construct
        LagrangianEqn(LagrangianEqn<Type>&& eqn);

        //- Construct from tmp
        LagrangianEqn(const tmp<LagrangianEqn<Type>>& tEqn);


    //- Destructor
    ~LagrangianEqn();


    // Member Functions

        //- Construct from another equation, with empty coefficients
        static tmp<LagrangianEqn<Type>> NewEmpty(const LagrangianEqn<Type>&);

        //- Return the field
        tmp<LagrangianSubSubField<Type>> psi() const;

        //- Return the field IO
        const regIOobject& psiIo() const;

        //- Return the field name
        const word& psiName() const;

        //- Return whether the given field is that of the equation
        template<template<class> class PrimitiveField>
        bool isPsi(const LagrangianSubField<Type, PrimitiveField>& psi) const;

        //- Check the operation with another equation
        void op(const LagrangianEqn<Type>& other);

        //- Set the previous field to that stored by another tmp equation
        void setPrevious(const tmp<LagrangianEqn<Type>>& tOther);

        //- Return the combined time and non-time explicit coefficient
        tmp<LagrangianCoeff<Type, false>> allSu() const;

        //- Return the combined time and non-time implicit coefficient
        tmp<LagrangianSp<Type>> allSp() const;

        //- Return the combined time and non-time explicit diagonal
        //  coefficient; i.e., with the off-diagonal parts of any tensor Sp
        //  coefficients converted into explicit values
        tmp<LagrangianCoeff<Type, false>> allDiagonalSu() const;

        //- Return the combined time and non-time implicit diagonal
        //  coefficient; i.e., with the off-diagonal parts of any tensor Sp
        //  coefficients removed
        tmp<LagrangianCoeff<scalar, true>> allDiagonalSp() const;

        //- Solve
        void solve(const bool final);


    // Member Operators

        //- Addition assignment
        void operator+=(const LagrangianEqn<Type>& other);

        //- Addition assignment
        void operator+=(const tmp<LagrangianEqn<Type>>& tOther);

        //- Subtraction assignment
        void operator-=(const LagrangianEqn<Type>& other);

        //- Subtraction assignment
        void operator-=(const tmp<LagrangianEqn<Type>>& tOther);
};


//- Negation
template<class Type>
tmp<Foam::LagrangianEqn<Type>> operator-(const LagrangianEqn<Type>& eqn);

//- Negation
template<class Type>
tmp<LagrangianEqn<Type>> operator-(const tmp<LagrangianEqn<Type>>& tEqn);

#define LAGRANGIAN_EQN_EQN_OPERATOR(Op)                                        \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const LagrangianEqn<Type>& b                                           \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tA,                                    \
        const LagrangianEqn<Type>& b                                           \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tA,                                    \
        const tmp<LagrangianEqn<Type>>& tB                                     \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianEqn<Type>& a,                                          \
        const tmp<LagrangianEqn<Type>>& tB                                     \
    );

//- Addition
LAGRANGIAN_EQN_EQN_OPERATOR(+)

//- Subtraction
LAGRANGIAN_EQN_EQN_OPERATOR(-)

//- Set-equal-to
LAGRANGIAN_EQN_EQN_OPERATOR(==)

#undef LAGRANGIAN_EQN_EQN_OPERATOR

#define LAGRANGIAN_EQN_FIELD_OPERATOR(Op, LagrangianSubField)                  \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const LagrangianSubField<Type>& field                                  \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const LagrangianSubField<Type>& field                                  \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const tmp<LagrangianSubField<Type>>& tField                            \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const tmp<LagrangianSubField<Type>>& tField                            \
    );

#define LAGRANGIAN_FIELD_EQN_OPERATOR(Op, LagrangianSubField)                  \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianSubField<Type>& field,                                 \
        const LagrangianEqn<Type>& eqn                                         \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianSubField<Type>& field,                                 \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianSubField<Type>>& tField,                           \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianSubField<Type>>& tField,                           \
        const LagrangianEqn<Type>& eqn                                         \
    );

//- Addition
LAGRANGIAN_EQN_FIELD_OPERATOR(+, LagrangianSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(+, LagrangianSubField)
LAGRANGIAN_EQN_FIELD_OPERATOR(+, LagrangianSubSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(+, LagrangianSubSubField)

//- Subtraction
LAGRANGIAN_EQN_FIELD_OPERATOR(-, LagrangianSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(-, LagrangianSubField)
LAGRANGIAN_EQN_FIELD_OPERATOR(-, LagrangianSubSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(-, LagrangianSubSubField)

//- Set-equal-to
LAGRANGIAN_EQN_FIELD_OPERATOR(==, LagrangianSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(==, LagrangianSubField)
LAGRANGIAN_EQN_FIELD_OPERATOR(==, LagrangianSubSubField)
LAGRANGIAN_FIELD_EQN_OPERATOR(==, LagrangianSubSubField)

#undef LAGRANGIAN_EQN_FIELD_OPERATOR
#undef LAGRANGIAN_FIELD_EQN_OPERATOR

#define LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(Op, LagrangianSubField)           \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const LagrangianSubField<scalar>& field                                \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const LagrangianSubField<scalar>& field                                \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianEqn<Type>>& tEqn,                                  \
        const tmp<LagrangianSubField<scalar>>& tField                          \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianEqn<Type>& eqn,                                        \
        const tmp<LagrangianSubField<scalar>>& tField                          \
    );

#define LAGRANGIAN_SCALAR_FIELD_EQN_OPERATOR(Op, LagrangianSubField)           \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianSubField<scalar>& field,                               \
        const LagrangianEqn<Type>& eqn                                         \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const LagrangianSubField<scalar>& field,                               \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianSubField<scalar>>& tField,                         \
        const tmp<LagrangianEqn<Type>>& tEqn                                   \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<LagrangianEqn<Type>> operator Op                                       \
    (                                                                          \
        const tmp<LagrangianSubField<scalar>>& tField,                         \
        const LagrangianEqn<Type>& eqn                                         \
    );

//- Multiplication
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(*, LagrangianSubField)
LAGRANGIAN_SCALAR_FIELD_EQN_OPERATOR(*, LagrangianSubField)
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(*, LagrangianSubSubField)
LAGRANGIAN_SCALAR_FIELD_EQN_OPERATOR(*, LagrangianSubSubField)

//- Division
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(/, LagrangianSubField)
LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR(/, LagrangianSubSubField)

#undef LAGRANGIAN_EQN_SCALAR_FIELD_OPERATOR
#undef LAGRANGIAN_SCALAR_FIELD_EQN_OPERATOR

/*---------------------------------------------------------------------------*\
                    Class LagrangianEqn::PsiRef Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
template<template<class> class PrimitiveField>
struct Foam::LagrangianEqn<Type>::PsiRef
{
    //- Primitive field types match. Return a valid reference.
    inline static const LagrangianSubField<Type, PrimitiveField>& ref
    (
        const LagrangianSubField<Type, PrimitiveField>& psi
    )
    {
        return psi;
    };

    //- Primitive field types do not match. Return a null reference.
    template<template<class> class OtherPrimitiveField>
    inline static const LagrangianSubField<Type, PrimitiveField>& ref
    (
        const LagrangianSubField<Type, OtherPrimitiveField>&
    )
    {
        return NullObjectRef<LagrangianSubField<Type, PrimitiveField>>();
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
