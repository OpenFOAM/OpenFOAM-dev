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
    Foam::distributionSizeGroupFvPatchScalarField

Description
    Boundary condition for size-group fractions that sets values corresponding
    to a given distribution

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        distributionSizeGroup;
        distribution
        {
            type        RosinRammler;
            min         1 [mm];
            max         15 [mm];
            d           14 [mm];
            n           2 [];
        }
        value       $internalField;
    }
    \endverbatim

SourceFiles
    distributionSizeGroupFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef distributionSizeGroupFvPatchScalarField_H
#define distributionSizeGroupFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class distributionSizeGroupFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class distributionSizeGroupFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Distribution
        const autoPtr<distribution> distribution_;

        //- Cached allocation coefficient value. Calculated on demand.
        mutable autoPtr<scalar> etaPtr_;


public:

    //- Runtime type information
    TypeName("distributionSizeGroup");


    // Constructors

        //- Construct from patch, internal field and dictionary
        distributionSizeGroupFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        distributionSizeGroupFvPatchScalarField
        (
            const distributionSizeGroupFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        distributionSizeGroupFvPatchScalarField
        (
            const distributionSizeGroupFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        distributionSizeGroupFvPatchScalarField
        (
            const distributionSizeGroupFvPatchScalarField&,
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
                new distributionSizeGroupFvPatchScalarField
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
