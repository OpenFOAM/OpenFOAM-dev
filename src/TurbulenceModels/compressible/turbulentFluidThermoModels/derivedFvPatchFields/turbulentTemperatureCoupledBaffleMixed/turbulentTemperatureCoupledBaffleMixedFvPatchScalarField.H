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
    Foam::compressible::
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles. Optional thin thermal layer resistances can be
    specified through thicknessLayers and kappaLayers entries.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - refGradient = zero gradient
      - refValue = neighbour value
      - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    where KDelta is heat-transfer coefficient K * deltaCoeffs

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivities per layer [W/m/K] | no |
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        kappaMethod     lookup;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::temperatureCoupledBase

SourceFiles
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_H
#define turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class turbulentTemperatureCoupledBaffleMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private Data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureCoupledBaffleMixed");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField&,
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
                new turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
