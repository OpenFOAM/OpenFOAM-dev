/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField

Description
    A thermal wall function for simulation of subcooled nucleate wall boiling
    with runtime selectable submodels for:
      - wall heat flux partitioning model
      - nucleation site density
      - bubble departure frequency
      - bubble departure diameter

    Implements a version of the well-known RPI wall boiling model
    (Kurul & Podowski, 1991). The model implementation based on implementation
    described in Peltola et al. (2019) and is similar to the model described by
    Peltola & Pättikangas (2012). The present implementation includes simplified
    support for presence of non-volatile components in addition to a single
    volatile component.

    References:
    \verbatim
        Kurul, N., & Podowski, M.Z. (1991).
        On the modeling of multidimensional effects in boiling channels.
        ANS. Proc. National Heat Transfer Con. Minneapolis, Minnesota, USA,
        1991.
        ISBN: 0-89448-162-1, pp. 30-40.
    \endverbatim

    \verbatim
        Peltola, J., Pättikangas, T., Bainbridge, W., Lehnigk, R., Schlegel, F.
        (2019).
        On Development and validation of subcooled nucleate boiling models for
        OpenFOAM Foundation Release.
        NURETH-18 Conference Proceedings, Portland, Oregon, United States, 2019.
    \endverbatim

    \verbatim
        Peltola, J., & Pättikangas, T.J.H. (2012).
        Development and validation of a boiling model for OpenFOAM multiphase
        solver.
        CFD4NRS-4 Conference Proceedings, Daejeon, Korea, 2012.
        paper 59.
    \endverbatim

Usage
    \table
        Property       | Description              | Required    | Default value
        phaseType      | 'vapour' or 'liquid'     | yes         |
        useLiquidTemperatureWallFunction | \\
        Use wall function to calculate liquid temperature? | no | yes
        tolerance      | solution tolerance       | no          | rootSmall
        Prt            | turbulent Prandtl number | no          | 0.85
    \endtable

        if phaseType 'vapour':
    \table
        partitioningModel |                       | yes         |
    \endtable

        if phaseType 'liquid':
    \table
        partitioningModel |                       | yes         |
        nucleationSiteModel |                     | yes         |
        departureDiameterModel |                  | yes         |
        departureFrequenctModel |                 | yes         |
        bubbleWaitingTimeRatio |                  | no          | 0.8
    \endtable

    NOTE: Runtime selectable submodels may require model specific entries

    Example usage:
    \verbatim
    hotWall
    {
        type            compressible::alphatWallBoilingWallFunction;
        phaseType       liquid;
        dmdt            uniform 0;
        Prt             0.85;
        partitioningModel
        {
            type        Lavieville;
            alphaCrit   0.2;
        }
        nucleationSiteModel
        {
            type        LemmertChawla;
        }
        departureDiameterModel
        {
            type        TolubinskiKostanchuk;
        }
        departureFrequencyModel
        {
            type        Cole;
        }
        value           uniform 0.01;
    \endverbatim

SourceFiles
    alphatWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressible_alphatWallBoilingWallFunctionFvPatchScalarField_H
#define compressible_alphatWallBoilingWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeWallFunctionBase.H"
#include "fixedValueFvPatchFields.H"
#include "partitioningModel.H"
#include "nucleationSiteModel.H"
#include "departureDiameterModel.H"
#include "departureFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
       Class alphatWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatWallBoilingWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField,
    public alphatPhaseChangeWallFunctionBase
{
public:

    // Data types

        //- Enumeration listing the possible operational modes
        enum phaseType
        {
            vapourPhase,
            vaporPhase, // <-- backwards compatibility
            liquidPhase
        };

        //- Heat source type names
        static const NamedEnum<phaseType, 3> phaseTypeNames_;


private:

    // Private classes

        //- Struct to hold cached properties for one of the phases
        struct properties;

        //- Struct to hold cached properties for the liquid when it is boiling
        struct boilingLiquidProperties;


    // Private Data

        // Controls

            //- Heat source type
            phaseType phaseType_;

            //- Estimate liquid temperature using logarithmic wall function?
            Switch useLiquidTemperatureWallFunction_;

            //- Solution tolerance
            scalar tolerance_;


        // Parameters

            //- Turbulent Prandtl number
            scalar Prt_;

            //- Bubble waiting time ratio
            scalar tau_;


        // Sub-Models

            //- Heat flux partitioning model
            autoPtr<wallBoilingModels::partitioningModel>
                partitioningModel_;

            //- Nucleation site density model
            autoPtr<wallBoilingModels::nucleationSiteModel>
                nucleationSiteModel_;

            //- Bubble departure diameter model
            autoPtr<wallBoilingModels::departureDiameterModel>
                departureDiameterModel_;

            //- Bubble departure frequency model
            autoPtr<wallBoilingModels::departureFrequencyModel>
                departureFrequencyModel_;


        // State

            //- Wall liquid fraction
            scalarField wetFraction_;

            //- Departure diameter
            scalarField dDeparture_;

            //- Departure frequency
            scalarField fDeparture_;

            //- Nucleation site density
            scalarField nucleationSiteDensity_;

            //- Quenching surface heat flux
            scalarField qQuenching_;

            //- Evaporative surface heat flux
            scalarField qEvaporative_;

            //- Mass transfer rate
            scalarField dmdtf_;


        // Private Member Functions

            //- Calculate the boiling for the given wall temperature. Return
            //  the total sum of all heat fluxes. Set the properties passed by
            //  non-const reference. Used by the functions below.
            tmp<scalarField> calcBoiling
            (
                const boilingLiquidProperties& props,
                const scalarField& Tw,
                scalarField& dDeparture,
                scalarField& fDeparture,
                scalarField& nucleationSiteDensity,
                scalarField& qQuenching,
                scalarField& qEvaporative,
                scalarField& dmdtf
            ) const;

            //- Calculate the boiling for the given wall temperature. Return
            //  the total sum of all heat fluxes. Use this to solve the balance
            //  between the heat fluxes specified by the boiling models and the
            //  temperature boundary condition without changing the stored
            //  boiling state.
            tmp<scalarField> calcBoiling
            (
                const boilingLiquidProperties& props,
                const scalarField& Tw
            ) const;

            //- Calculate the boiling for the given wall temperature. Return
            //  the total sum of all heat fluxes. Also set the stored boiling
            //  state. Use this after solving with the final wall temperature
            //  to set the boiling state.
            tmp<scalarField> evaluateBoiling
            (
                const boilingLiquidProperties& props,
                const scalarField& Tw
            );

            //- Get the temperature patch field and the parameters associated
            //  with its boundary condition that are necessary for
            //  approximately evaluating the boundary condition's heat flux at
            //  a given wall temperature.
            const fvPatchScalarField& getTemperaturePatchField
            (
                const boilingLiquidProperties& props,
                scalarField& isFixed,
                scalarField& h,
                scalarField& hTaPlusQa
            ) const;


public:

    //- Runtime type information
    TypeName("compressible::alphatWallBoilingWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatWallBoilingWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
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
                new alphatWallBoilingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Return the wall liquid fraction field [-]
        const scalarField& wetFraction() const
        {
            return wetFraction_;
        }

        //- Return the departure diameter field [m]
        const scalarField& dDeparture() const
        {
            return dDeparture_;
        }

        //- Return the departure frequency field [Hz]
        const scalarField& fDeparture() const
        {
            return fDeparture_;
        }

        //- Return the nucleation site density field [1/m^2]
        const scalarField& nucleationSiteDensity() const
        {
            return nucleationSiteDensity_;
        }

        //- Return the quenching surface heat flux field [W/m^2]
        const scalarField& qQuenching() const
        {
            return qQuenching_;
        }

        //- Return the evaporative surface heat flux field [W/m^2]
        const scalarField& qEvaporative() const
        {
            return qEvaporative_;
        }

        //- Return the rate of phase change
        virtual const scalarField& dmdtf() const
        {
            return dmdtf_;
        }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


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
