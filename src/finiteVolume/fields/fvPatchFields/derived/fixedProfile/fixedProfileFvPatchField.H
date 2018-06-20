/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::fixedProfileFvPatchField

Description
    This boundary condition provides a fixed value profile condition.

Usage
    \table
        Property     | Description              | Required | Default value
        profile      | Profile Function1        | yes |
        direction    | Profile direction        | yes |
        origin       | Profile origin           | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedProfile;
        profile         csvFile;

        profileCoeffs
        {
            nHeaderLine         0;          // Number of header lines
            refColumn           0;          // Reference column index
            componentColumns    (1 2 3);    // Component column indices
            separator           ",";        // Optional (defaults to ",")
            mergeSeparators     no;         // Merge multiple separators
            file                "Uprofile.csv";  // name of csv data file
            outOfBounds         clamp;      // Optional out-of-bounds handling
            interpolationScheme linear;     // Optional interpolation scheme
        }
        direction        (0 1 0);
        origin           0;
    }
    \endverbatim

    Example setting a parabolic inlet profile for the PitzDaily case:
    \verbatim
    inlet
    {
        type            fixedProfile;

        profile         polynomial
        (
            ((1 0 0)        (0 0 0))
            ((-6200 0 0)    (2 0 0))
        );
        direction       (0 1 0);
        origin          0.0127;
    }
    \endverbatim

Note
    The profile entry is a Function1 type.  The example above gives the
    usage for supplying csv file.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    Foam::timeVaryingMappedFixedValueFvPatchField

SourceFiles
    fixedProfileFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedProfileFvPatchField_H
#define fixedProfileFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fixedProfileFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedProfileFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Profile data
        autoPtr<Function1<Type>> profile_;

        //- Profile direction
        vector dir_;

        //- Profile origin
        scalar origin_;


public:

    //- Runtime type information
    TypeName("fixedProfile");


    // Constructors

        //- Construct from patch and internal field
        fixedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        fixedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        fixedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedProfileFvPatchField
        //  onto a new patch
        fixedProfileFvPatchField
        (
            const fixedProfileFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedProfileFvPatchField
        (
            const fixedProfileFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fixedProfileFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedProfileFvPatchField
        (
            const fixedProfileFvPatchField<Type>&,
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
                new fixedProfileFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixedProfileFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
