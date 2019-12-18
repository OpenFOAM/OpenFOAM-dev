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
    Foam::gradientEnergyFvPatchScalarField

Description
    This boundary condition provides a gradient condition for energy. This is
    selected when the corresponding temperature condition is zeroGradient,
    fixedGradient or gradientEnergyCalculatedTemperature.

    If the temperature condition is fixedGradient or zeroGradient, the
    parameters of the energy condition are obtained from a linearisation of the
    corresponding parameters of the temperature condition.

    If the temperature condition is gradientEnergyCalculatedTemperature, then
    the temperature condition explicitly provides an energy gradient field, so
    no linearisation is needed.

See also
    Foam::zeroGradientFvPatchScalarField
    Foam::fixedGradientFvPatchScalarField
    Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField

SourceFiles
    gradientEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef gradientEnergyFvPatchScalarField_H
#define gradientEnergyFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class gradientEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class gradientEnergyFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("gradientEnergy");


    // Constructors

        //- Construct from patch and internal field
        gradientEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gradientEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given gradientEnergyFvPatchScalarField
        // onto a new patch
        gradientEnergyFvPatchScalarField
        (
            const gradientEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        gradientEnergyFvPatchScalarField
        (
            const gradientEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gradientEnergyFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        gradientEnergyFvPatchScalarField
        (
            const gradientEnergyFvPatchScalarField&,
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
                new gradientEnergyFvPatchScalarField(*this, iF)
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
