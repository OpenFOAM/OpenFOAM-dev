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
    Foam::LagrangianPatchField

Description
    Base class for Lagrangian boundary conditions

SourceFiles
    LagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianPatchField_H
#define LagrangianPatchField_H

#include "LagrangianPatch.H"
#include "LagrangianFieldsFwd.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class LagrangianMesh;
class LagrangianSubMesh;

template<class Type>
class LagrangianPatchField;

template<class Type>
class calculatedLagrangianPatchField;

// Forward declaration of friend functions and operators
template<class Type>
Ostream& operator<<(Ostream&, const LagrangianPatchField<Type>&);

/*---------------------------------------------------------------------------*\
                         Class LagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LagrangianPatchField
{
    // Private Data

        //- Reference to patch
        const LagrangianPatch& patch_;


protected:

    // Protected Data

        //- Reference to the internal IO object
        const regIOobject& internalIo_;

        //- Reference to the internal field
        const LagrangianInternalDynamicField<Type>& internalField_;

        //- Reference to the non-dynamic internal field
        const LagrangianInternalField<Type>& internalNonDynamicField_;


public:

    //- Declare friendship with patch fields of different types
    template<class Type2>
    friend class LagrangianPatchField;


    // Public Type Definitions

        //- The patch type
        typedef LagrangianPatch Patch;

        //- The calculated patch field type
        typedef calculatedLagrangianPatchField<Type> Calculated;


    //- Runtime type information
    TypeName("LagrangianPatchField");


    //- Debug switch to disallow the use of genericLagrangianPatchField
    static int disallowGenericLagrangianPatchField;


    // Declare run-time constructor selection tables

        //- Select given patch and internal field
        declareRunTimeSelectionTable
        (
            autoPtr,
            LagrangianPatchField,
            LagrangianPatch,
            (
                const LagrangianPatch& p,
                const regIOobject& iIo
            ),
            (p, iIo)
        );

        //- Select given patch, internal field and dictionary
        declareRunTimeSelectionTable
        (
            autoPtr,
            LagrangianPatchField,
            dictionary,
            (
                const LagrangianPatch& p,
                const regIOobject& iIo,
                const dictionary& dict
            ),
            (p, iIo, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        LagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        LagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        LagrangianPatchField(const LagrangianPatchField<Type>&);

        //- Construct as copy setting the internal field reference
        LagrangianPatchField
        (
            const LagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct as copy setting the patch and the internal field reference
        LagrangianPatchField
        (
            const LagrangianPatchField<Type>&,
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const = 0;

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianPatchField<Type>> clone
        (
            const regIOobject&
        ) const = 0;

        //- Construct and return a clone onto a new patch
        //  and setting the internal field reference
        virtual autoPtr<LagrangianPatchField<Type>> clone
        (
            const LagrangianPatch&,
            const regIOobject&
        ) const
        {
            NotImplemented;
            return autoPtr<LagrangianPatchField<Type>>(nullptr);
        }


    // Selectors

        //- Return a pointer to a new LagrangianPatchField with a given type
        //  name, corresponding to a LagrangianPatch and internal field
        static autoPtr<LagrangianPatchField<Type>> New
        (
            const word& patchFieldType,
            const word& actualPatchType,
            const LagrangianPatch& p,
            const regIOobject&
        );

        //- Return a pointer to a new LagrangianPatchField with a given type
        //  name, corresponding to a LagrangianPatch and internal field
        static autoPtr<LagrangianPatchField<Type>> New
        (
            const word&,
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Return a pointer to a new LagrangianPatchField corresponding to a
        //  LagrangianPatch and internal field, with type and parameters
        //  specified in a dictionary
        static autoPtr<LagrangianPatchField<Type>> New
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedPointPatchField
        template<class Type2>
        static autoPtr<LagrangianPatchField<Type>>
        NewCalculatedType(const LagrangianPatchField<Type2>&);


    //- Destructor
    virtual ~LagrangianPatchField();


    // Member Functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const LagrangianPatch& patch() const;

            //- Return internal dimensions reference
            const dimensionSet& internalDimensions() const;

            //- Return internal field reference
            const LagrangianInternalDynamicField<Type>& internalField() const;

            //- Generate an internal sub field corresponding to this patch
            tmp<LagrangianSubSubField<Type>> internalSubField() const;

            //- Generate a primitive sub field corresponding to this patch
            SubField<Type> primitiveSubField() const;

            //- Return the type of the calculated form of LagrangianPatchField
            static const word& calculatedType();


        //- Reset to the given LagrangianPatchField
        virtual void reset(const LagrangianPatchField<Type>&);


        // Evaluation

            //- Initialise evaluation of the patch field
            virtual void initEvaluate
            (
                PstreamBuffers&,
                const LagrangianScalarInternalDynamicField& fraction
            );

            //- Evaluate the patch field
            virtual void evaluate
            (
                PstreamBuffers&,
                const LagrangianScalarInternalDynamicField& fraction
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Generate a sub-sub field corresponding to this patch
        tmp<LagrangianSubSubField<Type>> operator()() const;

        //- Define assignment and compound assignment operations
        #define MEMBER_OPERATOR(MemberOp, OtherType)                           \
            virtual void operator MemberOp                                     \
            (                                                                  \
                const LagrangianPatchField<OtherType>&                         \
            );                                                                 \
            virtual void operator MemberOp                                     \
            (                                                                  \
                const LagrangianSubField<OtherType>&                           \
            );                                                                 \
            virtual void operator MemberOp                                     \
            (                                                                  \
                const LagrangianSubSubField<OtherType>&                        \
            );                                                                 \
            virtual void operator MemberOp(const UList<OtherType>&);           \
            virtual void operator MemberOp(const OtherType&)
        MEMBER_OPERATOR(=, Type);
        MEMBER_OPERATOR(==, Type);
        MEMBER_OPERATOR(+=, Type);
        MEMBER_OPERATOR(-=, Type);
        MEMBER_OPERATOR(*=, scalar);
        MEMBER_OPERATOR(/=, scalar);
        #undef MEMBER_OPERATOR

        //- Ostream operator
        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const LagrangianPatchField<Type>&
        );
};


template<class Type>
const LagrangianPatchField<Type>& operator+
(
    const LagrangianPatchField<Type>& ptf,
    const Type&
)
{
    return ptf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianPatchFieldFunctions.H"

#ifdef NoRepository
    #include "LagrangianPatchField.C"
    #include "calculatedLagrangianPatchField.H"
#endif


#define addToLagrangianPatchFieldRunTimeSelection(                             \
    PatchTypeField, typePatchTypeField)                                        \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


#define addNullConstructableToLagrangianPatchFieldRunTimeSelection(            \
    PatchTypeField, typePatchTypeField)                                        \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        LagrangianPatch                                                        \
    );                                                                         \
    addToLagrangianPatchFieldRunTimeSelection                                  \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField                                                     \
    )


#define makeLagrangianPatchTypeField(                                          \
    PatchTypeField, typePatchTypeField)                                        \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToLagrangianPatchFieldRunTimeSelection                                  \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField                                                     \
    )


#define makeNullConstructableLagrangianPatchTypeField(                         \
    PatchTypeField, typePatchTypeField)                                        \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addNullConstructableToLagrangianPatchFieldRunTimeSelection                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField                                                     \
    )


#define makeTemplateLagrangianPatchTypeField(fieldType, type)                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, LagrangianPatch, CAPITALIZE(fieldType), Field),             \
        0                                                                      \
    );                                                                         \
    addToLagrangianPatchFieldRunTimeSelection                                  \
    (                                                                          \
        CAT3(LagrangianPatch, CAPITALIZE(fieldType), Field),                   \
        CAT4(type, LagrangianPatch, CAPITALIZE(fieldType), Field)              \
    )


#define makeNullConstructableTemplateLagrangianPatchTypeField(fieldType, type) \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, LagrangianPatch, CAPITALIZE(fieldType), Field),             \
        0                                                                      \
    );                                                                         \
    addNullConstructableToLagrangianPatchFieldRunTimeSelection                 \
    (                                                                          \
        CAT3(LagrangianPatch, CAPITALIZE(fieldType), Field),                   \
        CAT4(type, LagrangianPatch, CAPITALIZE(fieldType), Field)              \
    )


#define makeLagrangianPatchFields(type)                                        \
    FOR_ALL_FIELD_TYPES                                                        \
    (                                                                          \
        makeTemplateLagrangianPatchTypeField,                                  \
        type                                                                   \
    )


#define makeNullConstructableLagrangianPatchFields(type)                       \
    FOR_ALL_FIELD_TYPES                                                        \
    (                                                                          \
        makeNullConstructableTemplateLagrangianPatchTypeField,                 \
        type                                                                   \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
