/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::freestreamVelocityFvPatchVectorField

Description
    This boundary condition provides a free-stream condition for velocity.

    It is an inlet-outlet condition that uses the velocity orientation to
    continuously blend between fixed value for normal inlet and zero gradient
    for normal outlet flow.

Usage
    \table
        Property        | Description             | Required    | Default value
        freestreamValue | freestream velocity     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            freestreamVelocity;
        freestreamValue uniform (300 0 0);
    }
    \endverbatim

    Note:
        This condition is designed to operate with the freestreamPressure
        condition

See also
    Foam::mixedFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    freestreamVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef freestreamVelocityFvPatchVectorField_H
#define freestreamVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class freestreamVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class freestreamVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("freestreamVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        freestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given freestreamVelocityFvPatchVectorField onto
        //  a new patch
        freestreamVelocityFvPatchVectorField
        (
            const freestreamVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        freestreamVelocityFvPatchVectorField
        (
            const freestreamVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        freestreamVelocityFvPatchVectorField
        (
            const freestreamVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new freestreamVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

            const vectorField& freestreamValue() const
            {
                return refValue();
            }

            vectorField& freestreamValue()
            {
                return refValue();
            }


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
