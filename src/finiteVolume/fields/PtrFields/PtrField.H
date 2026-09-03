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

Class
    Foam::PtrField

Description
    Field of objects which are dynamically allocated and stored in a PtrList

SourceFiles
    PtrField.C

\*---------------------------------------------------------------------------*/

#ifndef PtrField_H
#define PtrField_H

#include "tmp.H"
#include "PtrList.H"
#include "expressionExpressionable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PtrField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PtrField
:
    public tmp<PtrField<Type>>::refCount,
    public PtrList<Type>
{

public:

    // Public Types

        //- Primitive type
        using pType = typename Type::pType;

        //- Primitive component type
        using pCmptType = typename Type::pCmptType;

        //- Component type
        using cmptType = typename Type::cmptType;

        //- Component Field type
        using cmptField = typename Type::cmptField;


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        PtrField();

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit PtrField(const label);

        //- Copy constructor
        PtrField(const PtrField<Type>&);

        //- Move constructor
        PtrField(PtrField<Type>&&);

        //- Construct as copy of tmp<PtrField>
        PtrField(const tmp<PtrField<Type>>&);

        //- Construct as copy of a PtrField<Type> with optional clone arguments
        template<class ... Args>
        explicit PtrField(const PtrField<Type>&, const Args& ... args);

        //- Construct as copy of a PtrList<Type> with optional clone arguments
        template<class ... Args>
        explicit PtrField(const PtrList<Type>&, const Args& ... args);

        //- Construct from an expression with optional New selector arguments
        template
        <
            class Expression,
            class ... Args,
            class = EnableIfExpressionable<Expression>
        >
        explicit PtrField(const Expression& e, const Args& ... args);

        //- Clone
        tmp<PtrField<Type>> clone() const;

        //- If the argument is a temporary return
        //  otherwise construct a copy by calling New for each element
        template<class ... Args>
        static tmp<PtrField<Type>> New
        (
            const tmp<PtrField<Type>>& tf,
            const Args& ... args
        );


    // Member Functions

        //- Negate this field
        void negate();

        //- Return a component field expression of the field
        inline auto component(const direction&) const;

        //- Replace a component field of the field
        void replace(const direction, const PtrField<cmptField>&);

        //- Replace a component field of the field from an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void replace(const direction, const Expression&);

        //- Replace a component field of the field
        void replace(const direction, const cmptType&);

        //- Return the field transpose expression
        //  (only defined for second rank tensors)
        inline auto T() const;

        //- Bound the field to the given lower bound
        void boundLower(const pType& lower);

        //- Bound the field to the given upper bound
        void boundUpper(const pType& upper);

        //- Bound the field within the given lower and upper bounds
        void bound(const pType& lower, const pType& upper);


    // Member Operators

        void operator=(const PtrField<Type>&);
        void operator=(PtrField<Type>&&);
        void operator=(const tmp<PtrField<Type>>&);

        //- Assign to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator=(const Expression& e);

        void operator=(const pType&);

        void operator+=(const PtrField<Type>&);
        void operator+=(const tmp<PtrField<Type>>&);

        //- Assign-add to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator+=(const Expression& e);

        void operator+=(const pType&);

        void operator-=(const PtrField<Type>&);
        void operator-=(const tmp<PtrField<Type>>&);

        //- Assign-subtract to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>

        void operator-=(const Expression& e);

        void operator-=(const pType&);

        void operator*=(const PtrField<cmptField>&);
        void operator*=(const tmp<PtrField<cmptField>>&);

        //- Assign-multiply to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator*=(const Expression& e);

        void operator*=(const pCmptType&);

        void operator/=(const PtrField<cmptField>&);
        void operator/=(const tmp<PtrField<cmptField>>&);

        //- Assign-divide to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator/=(const Expression& e);

        void operator/=(const pCmptType&);
};


template<class Type>
Ostream& operator<<(Ostream&, const PtrField<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrFieldExpression.H"

#ifdef NoRepository
    #include "PtrField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
