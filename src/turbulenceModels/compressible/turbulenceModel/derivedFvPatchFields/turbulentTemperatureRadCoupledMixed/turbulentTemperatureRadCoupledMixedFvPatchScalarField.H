/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::
    compressible::
    turbulentTemperatureRadCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer
    to be used for in multiregion cases. Optional thin thermal layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.

    The thermal conductivity, \c kappa, can either be retrieved from the
    mesh database using the \c lookup option, or from a \c solidThermo
    or \c fluidThermo thermophysical package.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        kappa        | thermal conductivity option | yes     |
        kappaName    | name of thermal conductivity field | no  |
        Tnbr         | name of the field    | no | T
        QrNbr      | name of the radiative flux in the nbr region | no | none
        Qr         | name of the radiative flux in this region | no | none
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivites per layer [W/m/K] | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            compressible::turbulentTemperatureRadCoupledMixed;
        Tnbr            T;
        kappa           lookup;
        KappaName       kappa;
        QrNbr           Qr; // or none. Name of Qr field on neighbour region
        Qr              Qr; // or none. Name of Qr field on local region
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4)
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    Note: kappa : heat conduction at patch. Gets supplied how to
            lookup/calculate
     kappa:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'fluidThermo' : use fluidThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use solidThermo kappa()
    - 'directionalSolidThermo' directionalKappa()

SourceFiles
    turbulentTemperatureRadCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureRadCoupledMixedFvPatchScalarField_H
#define turbulentTemperatureRadCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class turbulentTemperatureRadCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureRadCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private data

        //- Name of field on the neighbour region
        const word TnbrName_;

         //- Name of the radiative heat flux in the neighbout region
        const word QrNbrName_;

        //- Name of the radiative heat flux in local region
        const word QrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureRadCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const
            turbulentTemperatureRadCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureRadCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const turbulentTemperatureRadCoupledMixedFvPatchScalarField&,
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
                new turbulentTemperatureRadCoupledMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
