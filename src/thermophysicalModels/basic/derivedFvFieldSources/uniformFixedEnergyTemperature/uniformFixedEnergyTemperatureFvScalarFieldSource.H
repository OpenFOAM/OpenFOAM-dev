/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::uniformFixedEnergyTemperatureFvScalarFieldSource

Description
    This source condition is applied to the temperature field, but provides a
    uniform fixed energy into the energy equation

Usage
    \table
        Property       | Description                | Required | Default value
        uniformHe      | uniform energy value       | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformFixedEnergyTemperature;
        uniformHe       3700000;
    }
    \endverbatim

SourceFiles
    uniformFixedEnergyTemperatureFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedEnergyTemperatureFvScalarFieldSource_H
#define uniformFixedEnergyTemperatureFvScalarFieldSource_H

#include "energyCalculatedTemperatureFvScalarFieldSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class uniformFixedEnergyTemperatureFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class uniformFixedEnergyTemperatureFvScalarFieldSource
:
    public energyCalculatedTemperatureFvScalarFieldSource
{
private:

    // Private Data

        //- Uniform value
        autoPtr<Function1<scalar>> uniformHe_;


public:

    //- Runtime type information
    TypeName("uniformFixedEnergyTemperature");


    // Constructors

        //- Construct from internal field and dictionary
        uniformFixedEnergyTemperatureFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        uniformFixedEnergyTemperatureFvScalarFieldSource
        (
            const uniformFixedEnergyTemperatureFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new uniformFixedEnergyTemperatureFvScalarFieldSource(*this, iF)
            );
        }


    //- Destructor
    virtual ~uniformFixedEnergyTemperatureFvScalarFieldSource();


    // Member Functions

        //- Return the source energy value
        virtual tmp<scalarField> sourceHeValue(const fvSource&) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
