/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::Field

Description
    Generic templated field type.

SourceFiles
    FieldFunctions.H
    FieldFunctionsM.H
    FieldM.H
    Field.C
    FieldFunctions.C
    FieldFunctionsM.C

\*---------------------------------------------------------------------------*/

#ifndef Field_H
#define Field_H

#include "tmp.H"
#include "direction.H"
#include "VectorSpace.H"
#include "scalarList.H"
#include "labelList.H"
#include "expressionExpressionable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class Field;

//- Pre-declare related SubField type
template<class Type>
class SubField;

template<class Type>
void writeEntry(Ostream& os, const Field<Type>&);

class unitSet;
class dictionary;

//- Field types
template<class Type>
struct FieldTraits
{
    //- Primitive element type
    using pType = Type;

    //- Primitive element component type
    using pCmptType = typename pTraits<Type>::cmptType;

    //- Element component type (pCmptType)
    using cmptType = pCmptType;
};

//- Field<Field> types
template<class Type>
struct FieldTraits<Field<Type>>
{
    //- Primitive element type
    using pType = Type;

    //- Primitive element component type
    using pCmptType = typename pTraits<Type>::cmptType;

    //- Element component type (Field of pCmptType)
    using cmptType = Field<pCmptType>;
};


/*---------------------------------------------------------------------------*\
                            Class Field Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Field
:
    public tmp<Field<Type>>::refCount,
    public List<Type>
{
    // Private Member Functions

        //- Return reference to copy of this field if mapF is the same field
        //  otherwise return mapF
        const UList<Type>& copySelf
        (
            const UList<Type>& mapF,
            tmp<Field<Type>>& tmapFcpy
        ) const;


public:

    // Public Types

        //- Primitive type
        using pType = typename FieldTraits<Type>::pType;

        //- Primitive component type
        using pCmptType = typename FieldTraits<Type>::pCmptType;

        //- Component type
        using cmptType = typename FieldTraits<Type>::cmptType;

        //- Component Field type
        using cmptField = Field<cmptType>;

        //- subField type
        using subField = SubField<Type>;


    // Static Data Members

        static const char* const typeName;


    // Static Member Functions

        //- Return a null field
        inline static const Field<Type>& null()
        {
            return NullObjectRef<Field<Type>>();
        }


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        Field();

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit Field(const label);

        //- Construct given size and initial value
        Field(const label, const Type&);

        //- Construct given size and initialised to zero
        Field(const label, const zero);

        //- Copy constructor of a UList\<Type\>
        explicit Field(const UList<Type>&);

        //- Mover constructor transferring the List contents
        explicit Field(List<Type>&&);

        //- Copy constructor of a UIndirectList\<Type\>
        explicit Field(const UIndirectList<Type>&);

        //- Copy constructor
        Field(const Field<Type>&);

        //- Copy constructor or reuse as specified.
        Field(Field<Type>&, bool reuse);

        //- Move constructor transferring the Field contents
        Field(Field<Type>&&);

        //- Copy constructor of tmp<Field>
        Field(const tmp<Field<Type>>&);

        //- Construct by 1 to 1 mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing
        );

        //- Construct by 1 to 1 mapping from the given tmp field
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by interpolative mapping from the given tmp field
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct from Istream
        Field(Istream&);

        //- Construct from a dictionary entry
        Field(const word& keyword, const dictionary&, const label size);

        //- Construct from a dictionary entry with unit conversion
        Field
        (
            const word& keyword,
            const unitSet&,
            const dictionary&,
            const label size
        );

        //- Construct from an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        Field(const Expression& e);

        //- Clone
        tmp<Field<Type>> clone() const;

        //- Return a pointer to a new calculatedFvPatchFieldField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<Field<Type>> NewCalculatedType(const Field<Type2>& f)
        {
            return tmp<Field<Type>>(new Field<Type>(f.size()));
        }


    // Member Functions

        //- 1 to 1 map from the given field
        void map
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing
        );

        //- 1 to 1 map from the given tmp field
        void map
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing
        );

        //- Interpolative map from the given field
        void map
        (
            const UList<Type>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Interpolative map from the given tmp field
        void map
        (
            const tmp<Field<Type>>& tmapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- 1 to 1 reverse-map from the given field
        void rmap
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing
        );

        //- 1 to 1 reverse-map from the given tmp field
        void rmap
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing
        );

        //- Interpolative reverse map from the given field
        void rmap
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing,
            const UList<scalar>& weights
        );

        //- Interpolative reverse map from the given tmp field
        void rmap
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing,
            const UList<scalar>& weights
        );

        //- Reset the field values to the given field
        //  Equivalent to operator=
        void reset(const Field<Type>&);

        //- Reset the field values to the given field
        //  Equivalent to operator=
        void reset(const UList<Type>&);

        //- Negate this field
        void negate();

        //- Return a component field expression of the field
        auto component(const direction& d) const;

        //- Replace a component field of the field
        void replace(const direction, const UList<cmptType>&);

        //- Replace a component field of the field
        void replace(const direction, const tmp<Field<cmptType>>&);

        //- Replace a component field of the field
        void replace(const direction, const cmptType&);

        template<class VSForm>
        VSForm block(const label start) const;

        //- Return the field transpose expression
        //  (only defined for second rank tensors)
        auto T() const;

        //- Append an element at the end of the Field
        void append(const Type& t);

        //- Append a Field at the end of the Field
        void append(const Field<Type>& f);

        //- Append a Field at the end of the Field
        void append(const tmp<Field<Type>>& tf);

        //- Append a Field evaluated from an expression at the end of this Field
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void append(const Expression& e);


    // Member Operators

        void operator=(const Field<Type>&);
        void operator=(Field<Type>&&);
        void operator=(const UList<Type>&);
        void operator=(List<Type>&&);
        void operator=(const tmp<Field<Type>>&);
        void operator=(const pType&);
        void operator=(const zero);

        template<class Form, class Cmpt, direction nCmpt>
        void operator=(const VectorSpace<Form,Cmpt,nCmpt>&);

        //- Assign to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator=(const Expression& e);

        void operator+=(const UList<Type>&);
        void operator+=(const Field<Type>&);
        void operator+=(const tmp<Field<Type>>&);
        void operator+=(const pType&);

        //- Assign-add to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator+=(const Expression& e);

        void operator-=(const UList<Type>&);
        void operator-=(const Field<Type>&);
        void operator-=(const tmp<Field<Type>>&);
        void operator-=(const pType&);

        //- Assign-subtract to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator-=(const Expression& e);

        void operator*=(const UList<cmptType>&);
        void operator*=(const Field<cmptType>&);
        void operator*=(const tmp<Field<cmptType>>&);
        void operator*=(const pCmptType&);

        //- Assign-multiply to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator*=(const Expression& e);

        void operator/=(const UList<cmptType>&);
        void operator/=(const Field<cmptType>&);
        void operator/=(const tmp<Field<cmptType>>&);
        void operator/=(const pCmptType&);

        //- Assign-divide to an expression
        template<class Expression, class = EnableIfExpressionable<Expression>>
        void operator/=(const Expression& e);
};


template<class Type>
Ostream& operator<<(Ostream&, const Field<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFunctions.H"
#include "FieldExpression.H"
#include "scalarField.H"

#ifdef NoRepository
    #include "Field.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
