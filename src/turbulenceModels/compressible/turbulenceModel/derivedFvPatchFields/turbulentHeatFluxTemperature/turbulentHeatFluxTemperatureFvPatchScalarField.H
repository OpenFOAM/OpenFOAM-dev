/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::turbulentHeatFluxTemperatureFvPatchScalarField

Description
    Fixed heat boundary condition to specify temperature gradient. Input
    heat source either specified in terms of an absolute power [W], or as a
    flux [W/m2].

    Example usage:

        hotWall
        {
            type            compressible::turbulentHeatFluxTemperature;
            heatSource      flux;        // power [W]; flux [W/m2]
            q               uniform 10;  // heat power or flux
            kappa           fluidThermo; // calculate kappa=alphaEff*thermo.Cp
            Qr              none;        // name of the radiative flux
            value           uniform 300; // initial temperature value
        }


SourceFiles
    turbulentHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentHeatFluxTemperatureFvPatchScalarFields_H
#define turbulentHeatFluxTemperatureFvPatchScalarFields_H

#include "fixedGradientFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
      Class turbulentHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Data types

        //- Enumeration listing the possible hest source input modes
        enum heatSourceType
        {
            hsPower,
            hsFlux
        };


private:

    // Private data

        //- Heat source type names
        static const NamedEnum<heatSourceType, 2> heatSourceTypeNames_;

        //- Heat source type
        heatSourceType heatSource_;

        //- Heat power [W] or flux [W/m2]
        scalarField q_;

        //- Name of radiative in flux field
        word QrName_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentHeatFluxTemperatureFvPatchScalarField onto
        //  a new patch
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentHeatFluxTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentHeatFluxTemperatureFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureFvPatchScalarField&,
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
                new turbulentHeatFluxTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
