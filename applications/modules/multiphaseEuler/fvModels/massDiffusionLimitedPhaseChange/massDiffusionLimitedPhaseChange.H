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
    Foam::fv::massDiffusionLimitedPhaseChange

Description
    Model for mass-diffusion rate limited phase change between two phases.

    One of the two phases is considered to be the limiting phase with respect
    to the rate of mass-diffusion of a given specie. The transferring specie's
    composition on the limiting phase's side of the interface is given by a
    run-time selectable interface composition model. A diffusive mass transfer
    model then provides a coefficient which when multiplied by the difference
    in specie concentration between the interface and the bulk of the limiting
    phase gives the rate of mass transfer.

    All the specie transfers are combined and the associated latent heat is is
    equated to the rate of heat transfer from the two phases to the interface.
    This relation is solved for the interface state and the rate of phase
    change.

    This model requires at least one phase to be multi-component. A
    two-resistance heat transfer model must also be in operation between the
    two changing phases.

Usage
    Example usage:
    \verbatim
    phaseChange
    {
        type            massDiffusionLimitedPhaseChange;
        libs            ("libmultiphaseEulerFvModels.so");

        phases          (gas liquid);

        energySemiImplicit yes;

        interfaceComposition
        {
            gas_liquid_inThe_gas
            {
                type            saturated;
                species         (H2O);
                Le              1.0;
                pSat            ArdenBuck;
            }
        }

        diffusiveMassTransfer
        {
            blending        heatAndDiffusiveMassTransfer;

            gas_dispersedIn_liquid_inThe_gas
            {
                type            spherical;
            }

            liquid_dispersedIn_gas_inThe_gas
            {
                type            Frossling;
                Le              1.0;
            }
        }
    }
    \endverbatim

SourceFiles
    massDiffusionLimitedPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef massDiffusionLimitedPhaseChange_H
#define massDiffusionLimitedPhaseChange_H

#include "phaseSystem.H"
#include "phaseChange.H"
#include "interfaceCompositionModel.H"
#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace solvers
{
    class multiphaseEuler;
}

namespace fv
{

/*---------------------------------------------------------------------------*\
              Class massDiffusionLimitedPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class massDiffusionLimitedPhaseChange
:
    public phaseChange
{
private:

    // Private Data

        //- Solver
        const solvers::multiphaseEuler& solver_;

        //- Phase system
        const phaseSystem& fluid_;

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;

        //- The interface composition model
        autoPtr<sidedInterfaceCompositionModel> interfaceCompositionModel_;

        //- The diffusive mass transfer model
        autoPtr<blendedSidedDiffusiveMassTransferModel>
            diffusiveMassTransferModel_;

        //- The number of correctors applied to the solution of the surface
        //  temperature (i.e., the interfacial heat and mass transfer balance)
        label nIter_;

        //- The surface temperature
        volScalarField::Internal Ts_;

        //- Explicit coefficients for each side and each specie
        Pair<PtrList<volScalarField::Internal>> mDotSus_;

        //- Implicit coefficients for each side and each specie
        Pair<PtrList<volScalarField::Internal>> mDotSps_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //- The phase change rate
        mutable volScalarField::Internal mDot_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Combine the lists of species from the interface composition models
        wordList getSpecies() const;

        //- Correct the phase change rate
        void correctMDot() const;


public:

    //- Runtime type information
    TypeName("massDiffusionLimitedPhaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        massDiffusionLimitedPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Evaluation

            //- Return the temperature at which the phases are considered to be
            //  changing
            virtual tmp<DimensionedField<scalar, volMesh>> Tchange() const;

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the mass transfer rate of a specie
            virtual tmp<DimensionedField<scalar, volMesh>> mDot
            (
                const label mDoti
            ) const;

            //- Override the pressure equation to trigger correction of the
            //  phase change rate
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the species equations to linearise in the mass fraction
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& heOrYi,
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

#endif

// ************************************************************************* //
