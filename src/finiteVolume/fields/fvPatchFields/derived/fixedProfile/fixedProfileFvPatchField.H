/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
        Property     | Description       | Required | Default value
        profile      | Profile function  | yes      |
        direction    | Direction of the line along which the profile is \\
                               evaluated | yes      |
        origin       | Origin of the line along which the profile is \\
                               evaluated | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedProfile;
        profile         table;
        profileCoeffs
        {
            file                "UProfile";
            format              csv;
            nHeaderLine         0;
            columns             (0 (1 2 3));
            separator           ",";
            mergeSeparators     no;
            outOfBounds         clamp;
            interpolationScheme linear;
        }
        direction        (0 1 0);
        origin           0;
    }
    \endverbatim

    Example setting a parabolic inlet profile for the pitzDaily case:
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

See also
    Foam::Function1s

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
    // Private Data

        //- Profile function
        autoPtr<Function1<Type>> profile_;

        //- Origin of the line along which the profile is evaluated
        scalar origin_;

        //- Direction of the line along which the profile is evaluated
        vector direction_;


public:

    //- Runtime type information
    TypeName("fixedProfile");


    // Constructors

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
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fixedProfileFvPatchField
        (
            const fixedProfileFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
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


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchField<Type>::operator=;
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
