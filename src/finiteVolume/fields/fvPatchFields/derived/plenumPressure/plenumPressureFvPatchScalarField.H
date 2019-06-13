/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::plenumPressureFvPatchScalarField

Description
    This boundary condition provides a plenum pressure inlet condition. This
    condition creates a zero-dimensional model of an enclosed volume of gas
    upstream of the inlet. The pressure that the boundary condition exerts on
    the inlet boundary is dependent on the thermodynamic state of the upstream
    volume. The upstream plenum density and temperature are time-stepped along
    with the rest of the simulation, and momentum is neglected. The plenum is
    supplied with a user specified mass flow and temperature.

    The result is a boundary condition which blends between a pressure inlet
    condition condition and a fixed mass flow. The smaller the plenum
    volume, the quicker the pressure responds to a deviation from the supply
    mass flow, and the closer the model approximates a fixed mass flow. As
    the plenum size increases, the model becomes more similar to a specified
    pressure.

    The expansion from the plenum to the inlet boundary is controlled by an
    area ratio and a discharge coefficient. The area ratio can be used to
    represent further acceleration between a sub-grid blockage such as fins.
    The discharge coefficient represents a fractional deviation from an
    ideal expansion process.

    This condition is useful for simulating unsteady internal flow problems
    for which both a mass flow boundary is unrealistic, and a pressure
    boundary is susceptible to flow reversal. It was developed for use in
    simulating confined combustion.

    Reference:
    \verbatim
        Bainbridge, W. (2013).
        The Numerical Simulation of Oscillations in Gas Turbine Combustion
        Chambers,
        PhD Thesis,
        Chapter 4, Section 4.3.1.2, 77-80.
    \endverbatim

Usage
    \table
        Property        | Description                 | Required | Default value
        gamma           | ratio of specific heats     | yes      | none
        R               | specific gas constant       | yes      | none
        supplyMassFlowRate | flow rate into the plenum | yes     | none
        supplyTotalTemperature | temperature into the plenum | yes | none
        plenumVolume    | plenum volume               | yes      | none
        plenumDensity   | plenum density              | yes      | none
        plenumTemperature | plenum temperature        | yes      | none
        U               | velocity field name         | no       | U
        phi             | flux field name             | no       | phi
        rho             | inlet density               | no       | none
        inletAreaRatio  | inlet open fraction         | yes      | none
        inletDischargeCoefficient | inlet loss coefficient | yes | none
        timeScale       | relaxation time scale       | yes      | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            plenumPressure;
        gamma           1.4;
        R               287.04;
        supplyMassFlowRate 0.0001;
        supplyTotalTemperature 300;
        plenumVolume    0.000125;
        plenumDensity   1.1613;
        plenumTemperature 300;
        inletAreaRatio  1.0;
        inletDischargeCoefficient 0.8;
        timeScale       1e-4;
        value           uniform 1e5;
    }
    \endverbatim

SourceFiles
    plenumPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef plenumPressureFvPatchScalarField_H
#define plenumPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plenumPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class plenumPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Ratio of specific heats
        scalar gamma_;

        //- Specific gas constant
        scalar R_;

        //- Mass flow rate supplied to the plenum
        scalar supplyMassFlowRate_;

        //- Total temperature of the gas supplied to the plenum
        scalar supplyTotalTemperature_;

        //- The volume of the plenum
        scalar plenumVolume_;

        //- The mean density of the gas in the plenum
        scalar plenumDensity_;

        //- The old-time mean density of the gas in the plenum
        scalar plenumDensityOld_;

        //- The mean temperature of the gas in the plenum
        scalar plenumTemperature_;

        //- The mean old-time temperature of the gas in the plenum
        scalar plenumTemperatureOld_;

        //- The constant density used when phi is volumetric
        scalar rho_;

        //- Whether or not the constant density has been specified
        bool hasRho_;

        //- The ratio of open area to total area at the inlet
        //  Allows a grid or mesh to be represented
        scalar inletAreaRatio_;

        //- The discharge coefficient at the inlet
        scalar inletDischargeCoefficient_;

        //- The time scale over which changes in pressure are smoothed
        scalar timeScale_;

        //- The time index used for updating
        label timeIndex_;

        //- The name of the flux field
        word phiName_;

        //- The name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("plenumPressure");


    // Constructors

        //- Construct from patch and internal field
        plenumPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plenumPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given plenumPressureFvPatchScalarField
        //  onto a new patch
        plenumPressureFvPatchScalarField
        (
            const plenumPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        plenumPressureFvPatchScalarField
        (
            const plenumPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new plenumPressureFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        plenumPressureFvPatchScalarField
        (
            const plenumPressureFvPatchScalarField&,
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
                new plenumPressureFvPatchScalarField(*this, iF)
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
