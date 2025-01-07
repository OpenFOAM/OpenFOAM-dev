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
    Foam::LagrangianSp

Description
    Wrapper around LagrangianCoeff to specialise for the implicit coefficient.
    Trivial at present. The idea eventually is to make the vector implicit
    coefficient automatically build tensor coefficients as and when a model
    supplies tensor values.

\*---------------------------------------------------------------------------*/

#ifndef LagrangianSp_H
#define LagrangianSp_H

#include "LagrangianCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class LagrangianSp;

template<class Type>
class LagrangianEqn;

template<class Type>
tmp<LagrangianSubField<Type>> operator/
(
    const LagrangianCoeff<Type, false>& Su,
    const LagrangianSp<Type>& Sp
);

tmp<LagrangianSubField<vector>> operator/
(
    const LagrangianCoeff<vector, false>& Su,
    const LagrangianSp<vector>& Sp
);

/*---------------------------------------------------------------------------*\
                          Class LagrangianSp Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LagrangianSp
:
    public LagrangianCoeff<scalar, true>
{
public:

    // Constructors

        //- Inherit base class constructors
        using LagrangianCoeff<scalar, true>::LagrangianCoeff;

        //- Construct as copy
        LagrangianSp(const LagrangianSp<Type>& Sp);

        //- Construct as copy or re-use as specified
        LagrangianSp(LagrangianSp<Type>& Sp, const bool reuse);

        //- Move construct
        LagrangianSp(LagrangianSp<Type>&& Sp);


    // Member Functions

        //- Return the scalar diagonal coefficient
        tmp<LagrangianCoeff<scalar, true>> A() const;

        //- Return the scalar off-diagonal coefficients
        tmp<LagrangianCoeff<Type, false>> H() const;

        //- Return the equivalent explicit coefficient
        tmp<LagrangianCoeff<Type, false>> Su() const;

        //- Return the equivalent explicit coefficient. Pass the equation
        //  explicitly so this can be used in the equation destructor when the
        //  reference held by the coeff classes might no longer be valid.
        tmp<LagrangianCoeff<Type, false>> Su(const LagrangianEqn<Type>&) const;


    // Member Operators

        //- Addition assignment
        using LagrangianCoeff<scalar, true>::operator+=;

        //- Addition assignment
        void operator+=(const LagrangianSp<Type>& Sp);

        //- Subtraction assignment
        using LagrangianCoeff<scalar, true>::operator-=;

        //- Subtraction assignment
        void operator-=(const LagrangianSp<Type>& Sp);

        //- Multiply assignment
        using LagrangianCoeff<scalar, true>::operator*=;

        //- Division assignment
        using LagrangianCoeff<scalar, true>::operator/=;


    // Friend Operators

        //- Divide an Su by an Sp
        friend tmp<LagrangianSubField<Type>> operator/ <Type>
        (
            const LagrangianCoeff<Type, false>& Su,
            const LagrangianSp<Type>& Sp
        );
};


/*---------------------------------------------------------------------------*\
                     Class LagrangianSp<vector> Declaration
\*---------------------------------------------------------------------------*/

template<>
class LagrangianSp<vector>
:
    public tmp<LagrangianSp<vector>>::refCount
{
    // Private Member Data

        //- Scalar coefficient
        LagrangianCoeff<scalar, true> scalarCoeff_;

        //- Tensor coefficient
        LagrangianCoeff<tensor, true> tensorCoeff_;


    // Private Member Functions

        //- Promote the coefficient storage from scalar to tensor
        void makeTensor();


public:

    // Constructors

        //- Construct from a reference to the equation
        LagrangianSp(const LagrangianEqn<vector>& eqn);

        //- Construct as copy
        LagrangianSp(const LagrangianSp<vector>& coeff);

        //- Construct as copy or re-use as specified
        LagrangianSp(LagrangianSp<vector>& coeff, const bool reuse);

        //- Move construct
        LagrangianSp(LagrangianSp<vector>&& coeff);

        //- Construct from a field
        template<template<class> class PrimitiveField>
        LagrangianSp
        (
            const LagrangianEqnBase& eqn,
            const LagrangianSubField<scalar, PrimitiveField>&
        );

        //- Construct from a tmp field
        template<template<class> class PrimitiveField>
        LagrangianSp
        (
            const LagrangianEqnBase& eqn,
            const tmp<LagrangianSubField<scalar, PrimitiveField>>&
        );


    // Member Functions

        //- Access the equation
        const LagrangianEqnBase& eqn() const;

        //- Determine whether this coefficient has values or not
        bool valid() const;

        //- In-place negation
        void negate();

        //- Return the scalar diagonal coefficient
        tmp<LagrangianCoeff<scalar, true>> A() const;

        //- Return the scalar off-diagonal coefficients
        tmp<LagrangianCoeff<vector, false>> H() const;

        //- Return the equivalent explicit coefficient
        tmp<LagrangianCoeff<vector, false>> Su() const;

        //- Return the equivalent explicit coefficient. Pass the equation
        //  explicitly so this can be used in the equation destructor when the
        //  reference held by the coeff classes might no longer be valid.
        tmp<LagrangianCoeff<vector, false>> Su
        (
            const LagrangianEqn<vector>&
        ) const;


    // Member Operators

        //- Addition assignment
        template<template<class> class PrimitiveField>
        void operator+=(const LagrangianSubField<scalar, PrimitiveField>&);

        //- Addition assignment
        template<template<class> class PrimitiveField>
        void operator+=(const LagrangianSubField<tensor, PrimitiveField>&);

        //- Addition assignment
        template<template<class> class PrimitiveField>
        void operator+=(const tmp<LagrangianSubField<scalar, PrimitiveField>>&);

        //- Addition assignment
        template<template<class> class PrimitiveField>
        void operator+=(const tmp<LagrangianSubField<tensor, PrimitiveField>>&);

        //- Addition assignment
        void operator+=(const LagrangianCoeff<scalar, true>& coeff);

        //- Addition assignment
        void operator+=(const LagrangianCoeff<tensor, true>& coeff);

        //- Addition assignment
        void operator+=(const LagrangianSp<vector>& Sp);

        //- Addition assignment
        void operator+=(const dimensioned<scalar>&);

        //- Addition assignment
        void operator+=(const dimensioned<tensor>&);

        //- Addition assignment
        void operator+=(const zero);

        //- Subtraction assignment
        template<template<class> class PrimitiveField>
        void operator-=(const LagrangianSubField<scalar, PrimitiveField>&);

        //- Subtraction assignment
        template<template<class> class PrimitiveField>
        void operator-=(const LagrangianSubField<tensor, PrimitiveField>&);

        //- Subtraction assignment
        template<template<class> class PrimitiveField>
        void operator-=(const tmp<LagrangianSubField<scalar, PrimitiveField>>&);

        //- Subtraction assignment
        template<template<class> class PrimitiveField>
        void operator-=(const tmp<LagrangianSubField<tensor, PrimitiveField>>&);

        //- Addition assignment
        void operator-=(const LagrangianCoeff<scalar, true>& coeff);

        //- Addition assignment
        void operator-=(const LagrangianCoeff<tensor, true>& coeff);

        //- Subtraction assignment
        void operator-=(const LagrangianSp<vector>& Sp);

        //- Subtraction assignment
        void operator-=(const dimensioned<scalar>&);

        //- Subtraction assignment
        void operator-=(const dimensioned<tensor>&);

        //- Subtraction assignment
        void operator-=(const zero);

        //- Multiply assignment
        template<template<class> class PrimitiveField>
        void operator*=(const LagrangianSubField<scalar, PrimitiveField>&);

        //- Multiply assignment
        template<template<class> class PrimitiveField>
        void operator*=(const tmp<LagrangianSubField<scalar, PrimitiveField>>&);

        //- Multiply assignment
        void operator*=(const dimensioned<scalar>&);

        //- Division assignment
        template<template<class> class PrimitiveField>
        void operator/=(const LagrangianSubField<scalar, PrimitiveField>&);

        //- Division assignment
        template<template<class> class PrimitiveField>
        void operator/=(const tmp<LagrangianSubField<scalar, PrimitiveField>>&);

        //- Division assignment
        void operator/=(const dimensioned<scalar>& dt);


    // Friend Operators

        //- Divide an Su by an Sp
        friend tmp<LagrangianSubField<vector>> operator/
        (
            const LagrangianCoeff<vector, false>& Su,
            const LagrangianSp<vector>& Sp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianSp.C"
    #include "LagrangianVectorSpTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
