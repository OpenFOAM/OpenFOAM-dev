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
    Foam::externalTemperatureFvPatchScalarField

Description
    This boundary condition applies a heat flux condition to temperature
    on an external wall. Heat flux can be specified in the following ways:

      - fixed power: supply Q
      - fixed heat flux: supply q
      - fixed heat transfer coefficient: supply h and Ta

    where:
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
        h  | Heat transfer coefficient [W/m^2/K]
        Ta | Ambient temperature [K]
    \endvartable

    If more than one parameter is given then the heat fluxes are summed.

    If a heat transfer coefficient is given optional thin thermal layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.

    The patch thermal conductivity \c kappa is obtained from the region
    thermophysicalTransportModel so that this boundary condition can be applied
    directly to either fluid or solid regions.

    The ambient temperature Ta is specified as a Foam::Function1 of time but
    uniform is space.

Usage
    \table
    Property     | Description                 | Required | Default value
    Q            | Power [W]                   | no       | 0
    q            | Heat flux [W/m^2]           | no       | 0
    h            | Heat transfer coefficient [W/m^2/K] | no | 0
    Ta           | Ambient temperature [K]     | if h is given  |
    thicknessLayers | Layer thicknesses [m]    | no |
    kappaLayers  | Layer thermal conductivities [W/m/K] | no |
    relaxation   | Relaxation for the wall temperature | no | 1
    emissivity   | Surface emissivity for radiative flux to ambient | no | 0
    qr           | Name of the radiative field | no | none
    qrRelaxation | Relaxation factor for radiative field | no | 1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            externalTemperature;

        Ta              constant 300.0;
        h               uniform 10.0;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);

        value           $internalField;
    }
    \endverbatim

See also
    Foam::patchKappa
    Foam::mixedFvPatchScalarField

SourceFiles
    externalTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef externalTemperatureFvPatchScalarField_H
#define externalTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class externalTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class externalTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        // Heat power

            //- Do we have a heat power?
            bool haveQ_;

            //- Heat power [W]
            scalar Q_;


        // Heat flux

            //- Do we have a heat flux?
            bool haveq_;

            //- Heat flux [W/m^2]
            scalarField q_;


        // Heat transfer coefficient

            //- Do we have a heat transfer coefficient?
            bool haveh_;

            //- Heat transfer coefficient [W/m^2K]
            scalarField h_;

            //- Ambient temperature [K]
            autoPtr<Function1<scalar>> Ta_;

            //- Optional surface emissivity for radiative transfer to ambient
            scalar emissivity_;

            //- Thickness of layers
            scalarList thicknessLayers_;

            //- Conductivity of layers
            scalarList kappaLayers_;


        //- Relaxation for the wall temperature (thermal inertia)
        scalar relaxation_;


        // Radiation

            //- Name of the radiative heat flux
            const word qrName_;

            //- Relaxation for qr
            scalar qrRelaxation_;

            //- Cache qr for relaxation
            scalarField qrPrevious_;


public:

    //- Runtime type information
    TypeName("externalTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        externalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  externalTemperatureFvPatchScalarField
        //  onto a new patch
        externalTemperatureFvPatchScalarField
        (
            const externalTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        externalTemperatureFvPatchScalarField
        (
            const externalTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        externalTemperatureFvPatchScalarField
        (
            const externalTemperatureFvPatchScalarField&,
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
                new externalTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
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


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
