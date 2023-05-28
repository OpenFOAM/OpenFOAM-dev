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
    Foam::coupledTemperatureFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    with another region in a CHT case.  Optional thin wall material layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - refGradient = qs_/kappa
      - refValue = neighbour value
      - valueFraction = kappaByDeltaNbr/(kappaByDeltaNbr + kappaByDelta)

    where kappaByDelta is heat-transfer coefficient kappa*deltaCoeffs
    and qs is the optional source heat flux.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field       | no          | T
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivities per layer [W/m/K] | no |
        qs           | Optional source heat flux [W/m^2] | no | 0
        Qs           | Optional heat source [W]          | no | 0
        qrNbr        | name of the radiative flux in the nbr region | no | none
        qr           | name of the radiative flux in this region | no | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            coupledTemperature;
        thicknessLayers (0.1 0.2 0.3 0.4); // Optional wall layer thicknesses
        kappaLayers     (1 2 3 4); // Optional wall layer conductivities
        qs              uniform 100;    // Optional source heat flux [W/m^2]
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    The patch thermal conductivity \c kappa is obtained from the region
    thermophysicalTransportModel so that this boundary condition can be applied
    directly to either fluid or solid regions.

    Note that in order to provide an optional heat source either qs or Qs
    should be specified, not both.

SourceFiles
    coupledTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTemperatureFvPatchScalarField_H
#define coupledTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class coupledTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of the radiative heat flux in the neighbor region
        const word qrNbrName_;

        //- Name of the radiative heat flux in local region
        const word qrName_;

        //- Thickness of wall layers
        scalarList thicknessLayers_;

        //- Conductivity of wall layers
        scalarList kappaLayers_;

        //- Optional source heat flux [W/m^2]
        autoPtr<scalarField> qs_;

        //- Optional heat source [W]
        scalar Qs_;

        //- Optional wall material conductivity/thickness
        scalar wallKappaByDelta_;


protected:

    //- Get the patch kappa, kappa*Tc/delta and kappa/delta and also the
    //  heat-flux obtained from the sum heat-flux provided
    virtual void getThis
    (
        tmp<scalarField>& kappa,
        tmp<scalarField>& sumKappaTByDelta,
        tmp<scalarField>& sumKappaByDelta,
        scalarField& qTot,
        tmp<scalarField>& qByKappa
    ) const;

    //- Get the neighbour patch kappa*Tc/delta and kappa/delta
    //  and the heat-flux correction
    virtual void getNbr
    (
        tmp<scalarField>& sumKappaTByDeltaNbr,
        tmp<scalarField>& sumKappaByDeltaNbr,
        tmp<scalarField>& qNbr
    ) const;

    //- Get the neighbour patch wall temperature and heat-flux
    //  when wall thermal resistance is specified
    virtual void getNbr
    (
        tmp<scalarField>& TwNbr,
        tmp<scalarField>& qNbr
    ) const;

    //- Add field to result which may have not been previously set
    void add
    (
        tmp<scalarField>& result,
        const tmp<scalarField>& field
    ) const;


public:

    //- Runtime type information
    TypeName("coupledTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        coupledTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coupledTemperatureFvPatchScalarField onto a
        //  new patch
        coupledTemperatureFvPatchScalarField
        (
            const coupledTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coupledTemperatureFvPatchScalarField
        (
            const coupledTemperatureFvPatchScalarField&,
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
                new coupledTemperatureFvPatchScalarField
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
