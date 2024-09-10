/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::singleSizeGroupFvPatchScalarField

Description
    Boundary condition for size-group fractions that sets a value of one in a
    single group, identified by its index, and a value of zero in all other
    groups

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        singleSizeGroup;
        index       0;
        value       $internalField;
    }
    \endverbatim

SourceFiles
    singleSizeGroupFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef singleSizeGroupFvPatchScalarField_H
#define singleSizeGroupFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class singleSizeGroupFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class singleSizeGroupFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Size group index
        const scalar index_;


public:

    //- Runtime type information
    TypeName("singleSizeGroup");


    // Constructors

        //- Construct from patch, internal field and dictionary
        singleSizeGroupFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        singleSizeGroupFvPatchScalarField
        (
            const singleSizeGroupFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        singleSizeGroupFvPatchScalarField
        (
            const singleSizeGroupFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        singleSizeGroupFvPatchScalarField
        (
            const singleSizeGroupFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new singleSizeGroupFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
