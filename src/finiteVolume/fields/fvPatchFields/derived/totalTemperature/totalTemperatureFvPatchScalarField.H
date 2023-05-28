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
    Foam::totalTemperatureFvPatchScalarField

Description
    This boundary condition provides a total temperature condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | Velocity field name     | no          | U
        phi          | Flux field name         | no          | phi
        psi          | Compressibility field name | no       | psi
        gamma        | ratio of specific heats (Cp/Cv) | yes |
        T0           | reference temperature   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            totalTemperature;
        T0              uniform 300;
    }
    \endverbatim

SourceFiles
    totalTemperatureFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef totalTemperatureFvPatchScalarField_H
#define totalTemperatureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class totalTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class totalTemperatureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

        //- Reference temperature
        scalarField T0_;


public:

    //- Runtime type information
    TypeName("totalTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        totalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given totalTemperatureFvPatchScalarField
        //  onto a new patch
        totalTemperatureFvPatchScalarField
        (
            const totalTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        totalTemperatureFvPatchScalarField
        (
            const totalTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        totalTemperatureFvPatchScalarField
        (
            const totalTemperatureFvPatchScalarField&,
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
                new totalTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the reference temperature
            const scalarField& T0() const
            {
                return T0_;
            }

            //- Return reference to the temperature to allow adjustment
            scalarField& T0()
            {
                return T0_;
            }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchScalarField&,
                const fvPatchFieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


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
