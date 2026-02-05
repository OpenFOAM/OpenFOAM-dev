/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::fv::phaseSurfaceCondensation

Description
    Model for mass diffusion limited condensation on the surface of a
    third (solid) phase.

    This model functions very similarly to the wall condensation model (see that
    model for references). Compared to the wall condensation model, in this
    model the predicted condensation rate is assumed to be proportional to
    laminar diffusivity multiplied by a coefficient given by user selectable
    diffusive mass transfer model.

Usage
    Example usage:
    \verbatim
    phaseSurfaceCondensation
    {
        type            phaseSurfaceCondensation;
        libs            ("libmultiphaseEulerFvModels.so");

        // The (solid) phase
        phase           solid;

        // Note: Order is important. This model is one-way. It turns vapour
        // into liquid. The phases should be specified in this order.
        phases          (steam water);

        // The specie that is condensing
        specie  H2O;

        energySemiImplicit yes;
        specieSemiImplicit no;

        diffusiveMassTransfer
        {
            solid_dispersedIn_steam_inThe_steam
            {
                type            Frossling;
                Le              1;
            }
        }

        saturationPressure
        {
            type            ArdenBuck;
        }
    }
    \endverbatim

See also
    Foam::fv::wallCondensation
    Foam::fv::massDiffusionLimitedPhaseChange

SourceFiles
    phaseSurfaceCondensation.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSurfaceCondensation_H
#define phaseSurfaceCondensation_H

#include "phaseChange.H"
#include "phaseSystem.H"
#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace solvers
{
    class multiphaseEuler;
}

class saturationPressureModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class phaseSurfaceCondensation Declaration
\*---------------------------------------------------------------------------*/

class phaseSurfaceCondensation
:
    public phaseChange
{
    // Private Data

        //- Solver
        const solvers::multiphaseEuler& solver_;

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- Reference to the liquid phase
        const phaseModel& liquid_;

        //- Reference to the vapour phase
        const phaseModel& vapour_;

        //- Reference to the solid phase on which the boiling occurs
        const phaseModel& solid_;

        //- The diffusive mass transfer model
        autoPtr<blendedSidedDiffusiveMassTransferModel>
            diffusiveMassTransferModel_;

        //- The saturation curve
        autoPtr<saturationPressureModel> saturationModelPtr_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //-  Whether or not to linearise the specie source
        bool specieSemiImplicit_;

        //- Heat transfer rate
        mutable volScalarField::Internal q_;

        //- Mass transfer rate
        mutable volScalarField::Internal mDot_;

        //- The derivative of mass transfer rate
        mutable volScalarField::Internal mDotDy_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Correct the phase change rate
        void correctMDot() const;


public:

    //- Runtime type information
    TypeName("phaseSurfaceCondensation");


    // Constructors

        //- Construct from explicit source name and mesh
        //- Construct from components
        phaseSurfaceCondensation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Access

            //- Access the liquid phase
            inline const phaseModel& liquid() const;

            //- Access the vapour phase
            inline const phaseModel& vapour() const;

            //- Access the solid phase on which the boiling occurs
            inline const phaseModel& solid() const;


        // Evaluation

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the derivative of mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDotDy() const;

            //- Override the pressure equation to add the mass transfer rate
            //  linearised in the pressure
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the phase and vapours's energy equations to add
            //  additional latent heat and quenching heat transfer terms,
            //  handle optional linearisation of specie source
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

#include "phaseSurfaceCondensation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
