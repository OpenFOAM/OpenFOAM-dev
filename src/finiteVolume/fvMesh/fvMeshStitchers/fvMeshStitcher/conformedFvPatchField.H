/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::conformedFvPatchField

Description
    This vol field boundary condition holds data from both the original
    faces and any associated non-conformal faces, with the latter mapped to the
    conformal faces in the original patch. It is used during mesh change
    (between the un-stitch and stitch steps) to ensure that fields relating to
    both the original and the non-conformal patches are retained and mapped.

SourceFiles
    conformedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef conformedFvPatchField_H
#define conformedFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class conformedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class conformedFvPatchField
:
    public fvPatchField<Type>
{
    // Private Data

        //- The original patch field
        autoPtr<fvPatchField<Type>> origFieldPtr_;

        //- The non-conformal patch fields
        PtrList<fvPatchField<Type>> ncFieldPtrs_;

        //- The coverage for the non-conformal patch fields
        PtrList<scalarField> ncCoverages_;


    // Private Static Member Functions

        //- Map from the non-conformal patch index to the index stored within
        //  the conformed patch field on the original patch
        static labelList ncOrigNcField(const fvBoundaryMesh& fvbm);


    // Private Constructors

        //- Construct from components
        conformedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            autoPtr<fvPatchField<Type>>&& origFieldPtr,
            PtrList<fvPatchField<Type>>&& ncFieldPtrs,
            PtrList<scalarField>&& ncCoverages
        );


    // Private Member Functions

        //- Return the indices of the patches associated with the non-conformal
        //  patch fields
        labelList ncPatchIndices() const;



public:

    //- Runtime type information
    TypeName("conformed");


    // Static Member Functions

        //- Conform the given boundary field
        static void conform(typename VolField<Type>::Boundary& bF);

        //- Un-conform the given boundary field
        static void unconform(typename VolField<Type>::Boundary& bF);


    // Constructors

        //- Construct from patch and internal field
        conformedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        conformedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given conformedFvPatchField<Type>
        //  onto a new patch
        conformedFvPatchField
        (
            const conformedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        conformedFvPatchField
        (
            const conformedFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        conformedFvPatchField
        (
            const conformedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new conformedFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchField<Type>&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<Type>&);


        // Attributes

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return origFieldPtr_->fixesValue();
            }

            //- Return true if the value of the patch field
            //  is altered by assignment (the default)
            virtual bool assignable() const
            {
                return origFieldPtr_->assignable();
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return origFieldPtr_->coupled();
            }

            //- Return true if this overrides the underlying constraint type
            bool overridesConstraint() const
            {
                return origFieldPtr_->overridesConstraint();
            }


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type>> snGrad() const
            {
                return origFieldPtr_->snGrad();
            }

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {
                origFieldPtr_->initEvaluate(commsType);
            }

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            )
            {
                origFieldPtr_->evaluate(commsType);
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<Field<scalar>>& weights
            ) const
            {
                return origFieldPtr_->valueInternalCoeffs(weights);
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<Field<scalar>>& weights
            ) const
            {
                return origFieldPtr_->valueBoundaryCoeffs(weights);
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const
            {
                return origFieldPtr_->gradientInternalCoeffs();
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                return origFieldPtr_->gradientBoundaryCoeffs();
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "conformedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
