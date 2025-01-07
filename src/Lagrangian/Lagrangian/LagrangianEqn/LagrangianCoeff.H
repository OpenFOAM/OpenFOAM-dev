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
    Foam::LagrangianCoeff

Description
    Class to store a coefficient of a Lagrangian equation

SourceFiles
    LagrangianCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianCoeff_H
#define LagrangianCoeff_H

#include "LagrangianMesh.H"
#include "LagrangianSubFields.H"
#include "LagrangianEqnBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LagrangianCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type, bool Implicit>
class LagrangianCoeff
:
    public tmp<LagrangianCoeff<Type, Implicit>>::refCount
{
    // Private Data

        //- Reference to the equation
        const LagrangianEqnBase& eqn_;

        //- Coefficient
        autoPtr<LagrangianSubField<Type>> S_;


    // Private Member Functions

        //- Construct the coefficient field for given dimensions
        void initialise(const dimensionSet& dims);


public:

    // Constructors

        //- Construct from a reference to the equation
        LagrangianCoeff(const LagrangianEqnBase& eqn);

        //- Construct as copy
        LagrangianCoeff(const LagrangianCoeff<Type, Implicit>& coeff);

        //- Construct as copy or re-use as specified
        LagrangianCoeff
        (
            LagrangianCoeff<Type, Implicit>& coeff,
            const bool reuse
        );

        //- Move construct
        LagrangianCoeff(LagrangianCoeff<Type, Implicit>&& coeff);

        //- Construct from a field
        template<template<class> class PrimitiveField>
        LagrangianCoeff
        (
            const LagrangianEqnBase& eqn,
            const LagrangianSubField<Type, PrimitiveField>&
        );

        //- Construct from a tmp field
        template<template<class> class PrimitiveField>
        LagrangianCoeff
        (
            const LagrangianEqnBase& eqn,
            const tmp<LagrangianSubField<Type, PrimitiveField>>&
        );


    // Member Functions

        //- Access the equation
        const LagrangianEqnBase& eqn() const;

        //- Determine whether this coefficient has values or not
        bool valid() const;

        //- Access the coefficient
        const LagrangianSubField<Type>& S() const;

        //- In-place negation
        void negate();


    // Member Operators

        //- Addition assignment
        template<template<class> class PrimitiveField>
        void operator+=(const LagrangianSubField<Type, PrimitiveField>&);

        //- Addition assignment
        template<template<class> class PrimitiveField>
        void operator+=(const tmp<LagrangianSubField<Type, PrimitiveField>>&);

        //- Addition assignment
        void operator+=(const LagrangianCoeff<Type, Implicit>&);

        //- Addition assignment
        void operator+=(const dimensioned<Type>&);

        //- Addition assignment
        void operator+=(const zero);

        //- Subtraction assignment
        template<template<class> class PrimitiveField>
        void operator-=(const LagrangianSubField<Type, PrimitiveField>&);

        //- Subtraction assignment
        template<template<class> class PrimitiveField>
        void operator-=(const tmp<LagrangianSubField<Type, PrimitiveField>>&);

        //- Subtraction assignment
        void operator-=(const LagrangianCoeff<Type, Implicit>&);

        //- Subtraction assignment
        void operator-=(const dimensioned<Type>&);

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

        //- Multiply assignment
        void operator*=(const zero&);

        //- Division assignment
        template<template<class> class PrimitiveField>
        void operator/=(const LagrangianSubField<scalar, PrimitiveField>&);

        //- Division assignment
        template<template<class> class PrimitiveField>
        void operator/=(const tmp<LagrangianSubField<scalar, PrimitiveField>>&);

        //- Division assignment
        void operator/=(const dimensioned<scalar>& dt);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
