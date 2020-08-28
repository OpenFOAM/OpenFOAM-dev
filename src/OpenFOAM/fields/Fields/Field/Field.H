/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

template<class Type>
Ostream& operator<<(Ostream&, const Field<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const tmp<Field<Type>>&);

class dictionary;

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

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    //- Declare type of subField
    typedef SubField<Type> subField;


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

        //- Construct given size and initialed to zero
        Field(const label, const zero);

        //- Copy constructor of a UList\<Type\>
        explicit Field(const UList<Type>&);

        //- Mover constructor transferring the List contents
        explicit Field(List<Type>&&);

        //- Copy constructor of a UIndirectList\<Type\>
        explicit Field(const UIndirectList<Type>&);

        //- Copy constructor
        Field(const Field<Type>&);

        //- Copy constructor or re-use as specified.
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

        //- Negate this field
        void negate();

        //- Return a component field of the field
        tmp<Field<cmptType>> component(const direction) const;

        //- Replace a component field of the field
        void replace(const direction, const UList<cmptType>&);

        //- Replace a component field of the field
        void replace(const direction, const tmp<Field<cmptType>>&);

        //- Replace a component field of the field
        void replace(const direction, const cmptType&);

        template<class VSForm>
        VSForm block(const label start) const;

        //- Return the field transpose (only defined for second rank tensors)
        tmp<Field<Type>> T() const;


    // Member Operators

        void operator=(const Field<Type>&);
        void operator=(Field<Type>&&);
        void operator=(const UList<Type>&);
        void operator=(List<Type>&&);
        void operator=(const SubField<Type>&);
        void operator=(const tmp<Field<Type>>&);
        void operator=(const Type&);
        void operator=(const zero);

        template<class Form, class Cmpt, direction nCmpt>
        void operator=(const VectorSpace<Form,Cmpt,nCmpt>&);

        void operator+=(const UList<Type>&);
        void operator+=(const tmp<Field<Type>>&);

        void operator-=(const UList<Type>&);
        void operator-=(const tmp<Field<Type>>&);

        void operator*=(const UList<scalar>&);
        void operator*=(const tmp<Field<scalar>>&);

        void operator/=(const UList<scalar>&);
        void operator/=(const tmp<Field<scalar>>&);

        void operator+=(const Type&);
        void operator-=(const Type&);

        void operator*=(const scalar&);
        void operator/=(const scalar&);


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (Ostream&, const Field<Type>&);

        friend Ostream& operator<< <Type>
        (Ostream&, const tmp<Field<Type>>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFunctions.H"

#ifdef NoRepository
    #include "Field.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
