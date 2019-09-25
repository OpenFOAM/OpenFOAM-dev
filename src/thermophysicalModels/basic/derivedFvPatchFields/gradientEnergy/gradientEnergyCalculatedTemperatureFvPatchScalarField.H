/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::gradientEnergyCalculatedTemperatureFvPatchScalarField

Description
    Base class for temperature boundary conditions in which the parameters of
    the gradient energy condition can be set directly

SourceFiles
    gradientEnergyCalculatedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef gradientEnergyCalculatedTemperatureFvPatchScalarField_H
#define gradientEnergyCalculatedTemperatureFvPatchScalarField_H

#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class gradientEnergyCalculatedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class gradientEnergyCalculatedTemperatureFvPatchScalarField
:
    public calculatedFvPatchScalarField
{
    // Private Data

        //- Energy normal gradient field
        scalarField heGradient_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("gradientEnergyCalculatedTemperature");


    // Constructors

        //- Construct from patch and internal field
        gradientEnergyCalculatedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gradientEnergyCalculatedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        gradientEnergyCalculatedTemperatureFvPatchScalarField
        (
            const gradientEnergyCalculatedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        gradientEnergyCalculatedTemperatureFvPatchScalarField
        (
            const gradientEnergyCalculatedTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gradientEnergyCalculatedTemperatureFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        gradientEnergyCalculatedTemperatureFvPatchScalarField
        (
            const gradientEnergyCalculatedTemperatureFvPatchScalarField&,
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
                new gradientEnergyCalculatedTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return true if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Access the energy gradient
            inline scalarField& heGradient()
            {
                return heGradient_;
            }

            //- Access the energy gradient
            inline const scalarField& heGradient() const
            {
                return heGradient_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
