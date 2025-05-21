/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::fv::wallBoiling

Description
    Model for nucleate wall boiling between two phases on the surface of a
    number of wall patches.

    This model implements a version of the well-known RPI wall boiling model
    (Kurul & Podowski, 1991). The model is based on an implementation described
    in Peltola et al. (2019) and is similar to the model described by Peltola &
    Pättikangas (2012).

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
    Example usage:
    \verbatim
    wallBoiling
    {
        type            wallBoiling;
        libs            ("libmultiphaseEulerFvModels.so");

        // Note: Order is important. This model is one-way. It turns liquid
        // into vapour. The phases should be specified in this order.
        phases          (water steam);

        energySemiImplicit no;

        saturationTemperature
        {
            type            constant;
            value           372.76;
        }

        partitioningModel
        {
            type            Lavieville;
            alphaCrit       0.2;
        }

        nucleationSiteModel
        {
            type            LemmertChawla;
            Cn              1;
            NRef            30000000;
            deltaTRef       10;
        }

        departureDiameterModel
        {
            type            TolubinskiKostanchuk;
            dRef            0.00024;
            dMax            0.0014;
            dMin            1e-06;
        }

        departureFrequencyModel
        {
            type            KocamustafaogullariIshii;
            Cf              1.18;
        }
    }
    \endverbatim

    In addition to the above fvModel specification, wall patches on which
    boiling is to be calculated should have an alphatBoilingWallFunction
    boundary condition applied to the turbulent thermal diffusivity field.

See also
    Foam::alphatBoilingWallFunctionFvPatchScalarField

SourceFiles
    wallBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoiling_H
#define wallBoiling_H

#include "phaseChange.H"
#include "nucleation.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class saturationTemperatureModel;

namespace wallBoilingModels
{
    class partitioningModel;
    class nucleationSiteModel;
    class departureDiameterModel;
    class departureFrequencyModel;
}

class wallBoilingPhaseChangeRateFvPatchScalarField;

namespace fv
{

/*---------------------------------------------------------------------------*\
              Class wallBoiling Declaration
\*---------------------------------------------------------------------------*/

class wallBoiling
:
    public phaseChange,
    public nucleation
{
private:

    // Private classes

        //- Struct to hold properties that are held constant during the
        //  iterative solution procedure
        struct laggedProperties;


    // Private Data

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- Reference to the liquid phase
        const phaseModel& liquid_;

        //- Reference to the vapour phase
        const phaseModel& vapour_;

        //- Reference to the liquid turbulent thermal diffusivity
        const volScalarField& alphatLiquid_;

        //- Reference to the vapour turbulent thermal diffusivity
        const volScalarField& alphatVapour_;

        //- Reference to the field associated with the pressure equation
        const volScalarField& p_rgh_;

        //- Solution tolerance
        scalar tolerance_;

        //- Estimate liquid temperature using logarithmic wall function?
        bool liquidTemperatureWallFunction_;

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Bubble waiting time ratio
        scalar bubbleWaitingTimeRatio_;

        //- The saturation curve
        autoPtr<saturationTemperatureModel> saturationModelPtr_;

        //- Run-time selected heat flux partitioning model
        autoPtr<wallBoilingModels::partitioningModel> partitioningModel_;

        //- Run-time selected nucleation site density model
        autoPtr<wallBoilingModels::nucleationSiteModel> nucleationSiteModel_;

        //- Run-time selected bubble departure diameter model
        autoPtr<wallBoilingModels::departureDiameterModel>
            departureDiameterModel_;

        //- Run-time selected bubble departure frequency model
        autoPtr<wallBoilingModels::departureFrequencyModel>
            departureFrequencyModel_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //- The phase change rate
        mutable volScalarField mDot_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Return the boundary condition types for the phase change rate
        wordList mDotBoundaryTypes() const;

        //- Get the liquid temperature patch field and the parameters
        //  associated with its boundary condition that are necessary for
        //  approximately evaluating the boundary condition's heat flux at
        //  a given wall temperature.
        const fvPatchScalarField& getLiquidTemperaturePatchField
        (
            const laggedProperties& lagProps,
            scalarField& isFixed,
            scalarField& h,
            scalarField& hTaPlusQa
        ) const;

        //- Calculate the boiling for the given wall temperature. Return
        //  the total sum of all heat fluxes. Set the properties passed by
        //  non-const reference. Used by the functions below.
        tmp<scalarField> calcBoiling
        (
            const laggedProperties& lagProps,
            const scalarField& TLiquid,
            const scalarField& wetFraction,
            scalarField& dDeparture,
            scalarField& fDeparture,
            scalarField& nucleationSiteDensity,
            scalarField& qQuenching,
            scalarField& qEvaporative,
            scalarField& mDot
        ) const;

        //- Calculate the boiling for the given wall temperature. Return
        //  the total sum of all heat fluxes. Use this to solve the balance
        //  between the heat fluxes specified by the boiling models and the
        //  temperature boundary condition without changing the stored
        //  boiling state.
        tmp<scalarField> calcBoiling
        (
            const wallBoilingPhaseChangeRateFvPatchScalarField& mDot,
            const laggedProperties& lagProps,
            const scalarField& TLiquid
        ) const;

        //- Calculate the boiling for the given wall temperature. Return
        //  the total sum of all heat fluxes. Also set the stored boiling
        //  state. Use this after solving with the final wall temperature
        //  to set the boiling state.
        tmp<scalarField> evaluateBoiling
        (
            wallBoilingPhaseChangeRateFvPatchScalarField& mDot,
            const laggedProperties& lagProps,
            const scalarField& TLiquid
        ) const;

        //- Correct the phase change rate
        void correctMDot() const;


public:

    //- Runtime type information
    TypeName("wallBoiling");


    // Constructors

        //- Construct from explicit source name and mesh
        wallBoiling
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Access the liquid phase
            inline const phaseModel& liquid() const;

            //- Access the vapour phase
            inline const phaseModel& vapour() const;

            //- Access the liquid turbulent thermal diffusivity
            inline const volScalarField& alphatLiquid() const;

            //- Access the vapour turbulent thermal diffusivity
            inline const volScalarField& alphatVapour() const;

            //- Is the given patch boiling?
            bool isBoiling(const label patchi) const;


        // Evaluation

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;

            //- Return the diameter of nuclei
            virtual tmp<DimensionedField<scalar, volMesh>> d() const;

            //- Return the number rate at which nuclei are generated
            virtual tmp<DimensionedField<scalar, volMesh>> nDot() const;

            //- Return the nucleation time scale
            virtual tmp<DimensionedField<scalar, volMesh>> tau() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the mass transfer rate for the given patch
            const wallBoilingPhaseChangeRateFvPatchScalarField& mDotPf
            (
                const label patchi
            ) const;

            //- Return the mass transfer rate for the given patch
            wallBoilingPhaseChangeRateFvPatchScalarField& mDotPfRef
            (
                const label patchi
            ) const;

            //- Use phaseChange's source functions
            using phaseChange::addSup;

            //- Override the pressure equation to add the mass transfer rate
            //  linearised in the pressure
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;


        //- Correct the fvModel
        virtual void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallBoilingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
