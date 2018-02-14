/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fvPatchField

Description
    Abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised types with
    specific evaluation procedures, particularly with respect to specific
    fields.

SourceFiles
    fvPatchField.C
    fvPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchField_H
#define fvPatchField_H

#include "fvPatch.H"
#include "DimensionedField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class volMesh;


// Forward declaration of friend functions and operators

template<class Type>
class fvPatchField;

template<class Type>
class calculatedFvPatchField;

template<class Type>
class fvMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const fvPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                        Class fvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, volMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;

        //- Update index used so that manipulateMatrix is called only once
        //  during the construction of the matrix
        bool manipulatedMatrix_;

        //- Optional patch type, used to allow specified boundary conditions
        //  to be applied to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;


public:

    typedef fvPatch Patch;
    typedef calculatedFvPatchField<Type> Calculated;


    //- Runtime type information
    TypeName("fvPatchField");

    //- Debug switch to disallow the use of genericFvPatchField
    static int disallowGenericFvPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patchMapper,
            (
                const fvPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and value
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Type& value
        );

        //- Construct from patch and internal field and patch type
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const word& patchType
        );

        //- Construct from patch and internal field and patch field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given fvPatchField onto a new patch
        fvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvPatchField(const fvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>(new fvPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        fvPatchField
        (
            const fvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>(new fvPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvPatchField<Type>> New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static tmp<fvPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvPatchField mapped onto a new patch
        static tmp<fvPatchField<Type>> New
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        static tmp<fvPatchField<Type>> NewCalculatedType
        (
            const fvPatch&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<fvPatchField<Type>> NewCalculatedType
        (
            const fvPatchField<Type2>&
        );


    //- Destructor
    virtual ~fvPatchField<Type>()
    {}


    // Member functions

        // Attributes

            //- Return the type of the calculated for of fvPatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if the value of the patch field
            //  is altered by assignment (the default)
            virtual bool assignable() const
            {
                return true;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }


        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, volMesh>&
            internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Optional patch type
            const word& patchType() const
            {
                return patchType_;
            }

            //- Optional patch type
            word& patchType()
            {
                return patchType_;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }

            //- Return true if the matrix has already been manipulated
            bool manipulatedMatrix() const
            {
                return manipulatedMatrix_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type>> snGrad() const;

            //- Return patch-normal gradient for coupled-patches
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> snGrad
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            //  with a weight field (0..1). This weight field is usually
            //  provided as the amount of geometric overlap for 'duplicate'
            //  patches. Sets Updated to true
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Return internal field next to patch as patch field
            virtual void patchInternalField(Field<Type>&) const;

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField() const
            {
                NotImplemented;
                return *this;
            }

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this coupled patchField
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this coupled patchField
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
                return *this;
            }


            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<Type>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;

            //- Helper function to write the keyword and entry only if the
            //  values are not equal. The value is then output as value2
            template<class EntryType>
            void writeEntryIfDifferent
            (
                Ostream& os,
                const word& entryName,
                const EntryType& value1,
                const EntryType& value2
            ) const;


        // Check

            //- Check fvPatchField<Type> against given fvPatchField<Type>
            void check(const fvPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const fvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvPatchField.C"
    #include "calculatedFvPatchField.H"
#endif


#define makeFvPatchField(fvPatchTypeField)                                     \
                                                                               \
defineNamedTemplateTypeNameAndDebug(fvPatchTypeField, 0);                      \
template<>                                                                     \
int fvPatchTypeField::disallowGenericFvPatchField                              \
(                                                                              \
    debug::debugSwitch("disallowGenericFvPatchField", 0)                       \
);                                                                             \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patch);                  \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patchMapper);            \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, dictionary);


#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)    \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// Use with caution
#define addRemovableToPatchFieldRunTimeSelection\
(PatchTypeField, typePatchTypeField)                                           \
                                                                               \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// For non-templated patch fields
#define makePatchTypeField(PatchTypeField, typePatchTypeField)                 \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// For non-templated patch fields - use with caution
#define makeRemovablePatchTypeField(PatchTypeField, typePatchTypeField)        \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addRemovableToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// For templated patch fields
#define makeTemplatePatchTypeField(fieldType, type)                            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field),                     \
        0                                                                      \
    );                                                                         \
    addToPatchFieldRunTimeSelection                                            \
    (                                                                          \
        CAT3(fvPatch, CAPITALIZE(fieldType), Field),                           \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field)                      \
    )

#define makePatchFields(type)                                                  \
    FOR_ALL_FIELD_TYPES(makeTemplatePatchTypeField, type)

#define makePatchFieldTypeName(fieldType, type)                                \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field),                     \
        0                                                                      \
    );

#define makePatchFieldTypeNames(type)                                          \
    FOR_ALL_FIELD_TYPES(makePatchFieldTypeName, type)

#define makePatchTypeFieldTypedef(fieldType, type)                             \
    typedef type##FvPatchField<fieldType>                                      \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field);

#define makePatchTypeFieldTypedefs(type)                                       \
    FOR_ALL_FIELD_TYPES(makePatchTypeFieldTypedef, type)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
