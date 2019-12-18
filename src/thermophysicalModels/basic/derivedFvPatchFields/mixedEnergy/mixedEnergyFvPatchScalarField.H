/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::mixedEnergyFvPatchScalarField

Description
    This boundary condition provides a mixed condition for energy. This is
    selected when the corresponding temperature condition is mixed, or
    mixedEnergyCalculatedTemperature.

    If the temperature condition is mixed, the parameters of the energy
    condition are obtained from a linearisation of the corresponding parameters
    of the temperature condition.

    If the temperature condition is mixedEnergyCalculatedTemperature, then the
    temperature condition explicitly provides energy
    refValue/refGrad/valueFraction fields, so no linearisation is needed.

See also
    Foam::mixedFvPatchField
    Foam::mixedEnergyCalculatedTemperatureFvPatchField

SourceFiles
    mixedEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedEnergyFvPatchScalarField_H
#define mixedEnergyFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class mixedEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixedEnergyFvPatchScalarField
:
    public mixedFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("mixedEnergy");


    // Constructors

        //- Construct from patch and internal field
        mixedEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedEnergyFvPatchScalarField
        // onto a new patch
        mixedEnergyFvPatchScalarField
        (
            const mixedEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        mixedEnergyFvPatchScalarField
        (
            const mixedEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixedEnergyFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        mixedEnergyFvPatchScalarField
        (
            const mixedEnergyFvPatchScalarField&,
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
                new mixedEnergyFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
