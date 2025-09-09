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
    Foam::fv::heatTransferLimitedPhaseChange

Description
    Model for heat transfer rate limited phase change between two phases.

    The interface between the two phases is assumed to be at a saturated
    condition. This allows the temperature of the interface to be evaluated
    from a user-supplied saturation curve. This temperature then defines the
    heat flux being transferred to the interface from the surrounding fluid.
    The imbalance in the heat fluxes on either side of the interface is then
    divided by the latent heat of phase change in order to get the rate at
    which mass is being changed from one phase to the other.

    This model only supports pure phases. A two-resistance heat transfer model
    must also be in operation between the two changing phases.

Usage
    Example usage:
    \verbatim
    phaseChange
    {
        type            heatTransferLimitedPhaseChange;
        libs            ("libmultiphaseEulerFvModels.so");

        phases          (steam water);

        energySemiImplicit yes;
        pressureImplicit no;

        saturationTemperature
        {
            type            constant;
            value           372.76;
        }
    }
    \endverbatim

SourceFiles
    heatTransferLimitedPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferLimitedPhaseChange_H
#define heatTransferLimitedPhaseChange_H

#include "phaseChange.H"
#include "phaseSystem.H"
#include "saturationTemperatureModel.H"

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
              Class heatTransferLimitedPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class heatTransferLimitedPhaseChange
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

        //- The saturation curve
        autoPtr<saturationTemperatureModel> saturationModelPtr_;

        //- Should this phase change be linearised in the pressure equation?
        bool pressureImplicit_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //- The phase change rate
        mutable volScalarField::Internal mDot_;

        //- The derivative of the phase change rate w.r.t. the pressure. Only
        //  valid if pressure implicit is selected.
        mutable autoPtr<volScalarField::Internal> dmDotdpPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Correct the phase change rate
        void correctMDot() const;


public:

    //- Runtime type information
    TypeName("heatTransferLimitedPhaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        heatTransferLimitedPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Evaluation

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

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

#endif

// ************************************************************************* //
