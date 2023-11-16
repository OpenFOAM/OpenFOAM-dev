/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::conformedFvsPatchField

Description
    This surface field boundary condition holds data from both the original
    faces and any associated non-conformal faces, with the latter mapped to the
    conformal faces in the original patch. It is used during mesh change
    (between the un-stitch and stitch steps) to ensure that fields relating to
    both the original and the non-conformal patches are retained and mapped.

SourceFiles
    conformedFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef conformedFvsPatchField_H
#define conformedFvsPatchField_H

#include "fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class conformedFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class conformedFvsPatchField
:
    public fvsPatchField<Type>
{
    // Private Data

        //- The original patch field
        autoPtr<fvsPatchField<Type>> origFieldPtr_;

        //- The associated non-conformal patch field
        autoPtr<calculatedFvsPatchField<Type>> ncFieldPtr_;


    // Private Constructors

        //- Construct from components
        conformedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            autoPtr<fvsPatchField<Type>>&& origFieldPtr,
            autoPtr<calculatedFvsPatchField<Type>>&& ncFieldPtr
        );


public:

    //- Runtime type information
    TypeName("conformed");


    // Static Member Functions

        //- Conform the given boundary field
        static void conform(typename SurfaceField<Type>::Boundary& bF);

        //- Un-conform the given boundary field
        static void unconform(typename SurfaceField<Type>::Boundary& bF);


    // Constructors

        //- Construct from patch and internal field
        conformedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        conformedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given conformedFvsPatchField<Type>
        //  onto a new patch
        conformedFvsPatchField
        (
            const conformedFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        conformedFvsPatchField
        (
            const conformedFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        conformedFvsPatchField
        (
            const conformedFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new conformedFvsPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvsPatchField onto this fvsPatchField
            virtual void map(const fvsPatchField<Type>&, const fieldMapper&);

            //- Reset the fvsPatchField to the given fvsPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvsPatchField<Type>&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "conformedFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
