/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::fv::phaseSurfaceBoiling

Description
    Model for nucleate wall boiling on the surface of a third (solid) phase.

    This model functions very similarly to the wall boiling model (see that
    model for references). The same sub-models are used, with exactly the same
    specification syntax. The only difference is that the third phase must be
    additionally specified, and that a two-resistance heat transfer model must
    be in operation between the liquid phase and the third phase.

Usage
    Example usage:
    \verbatim
    phaseSurfaceBoiling
    {
        type            phaseSurfaceBoiling;
        libs            ("libmultiphaseEulerFvModels.so");

        phase           solid;

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
        }

        departureDiameterModel
        {
            type            TolubinskiKostanchuk;
        }

        departureFrequencyModel
        {
            type            KocamustafaogullariIshii;
            Cf              1.18;
        }
    }
    \endverbatim

See also
    Foam::fv::wallBoiling

SourceFiles
    wallBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSurfaceBoiling_H
#define phaseSurfaceBoiling_H

#include "phaseChange.H"
#include "nucleation.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace solvers
{
    class multiphaseEuler;
}

class saturationTemperatureModel;

namespace wallBoilingModels
{
    class partitioningModel;
    class nucleationSiteModel;
    class departureDiameterModel;
    class departureFrequencyModel;
}

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class phaseSurfaceBoiling Declaration
\*---------------------------------------------------------------------------*/

class phaseSurfaceBoiling
:
    public phaseChange,
    public nucleation
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

        //- The saturation curve
        autoPtr<saturationTemperatureModel> saturationModelPtr_;

        //- Run-time selected heat flux partitioning model
        autoPtr<Foam::wallBoilingModels::partitioningModel>
            partitioningModel_;

        //- Run-time selected nucleation site density model
        autoPtr<Foam::wallBoilingModels::nucleationSiteModel>
            nucleationSiteModel_;

        //- Run-time selected bubble departure diameter model
        autoPtr<Foam::wallBoilingModels::departureDiameterModel>
            departureDiameterModel_;

        //- Run-time selected bubble departure frequency model
        autoPtr<Foam::wallBoilingModels::departureFrequencyModel>
            departureFrequencyModel_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //- Surface liquid fraction
        mutable volScalarField::Internal wetFraction_;

        //- Departure diameter
        mutable volScalarField::Internal dDeparture_;

        //- Departure frequency
        mutable volScalarField::Internal fDeparture_;

        //- Nucleation site density
        mutable volScalarField::Internal nucleationSiteDensity_;

        //- Quenching heat transfer rate
        mutable volScalarField::Internal qQuenching_;

        //- Evaporative heat transfer rate
        mutable volScalarField::Internal qEvaporative_;

        //- Mass transfer rate
        mutable volScalarField::Internal mDot_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Correct the phase change rate
        void correctMDot() const;


public:

    //- Runtime type information
    TypeName("phaseSurfaceBoiling");


    // Constructors

        //- Construct from explicit source name and mesh
        //- Construct from components
        phaseSurfaceBoiling
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

            //- Return the diameter of nuclei
            virtual tmp<DimensionedField<scalar, volMesh>> d() const;

            //- Return the number rate at which nuclei are generated
            virtual tmp<DimensionedField<scalar, volMesh>> nDot() const;

            //- Return the nucleation time scale
            virtual tmp<DimensionedField<scalar, volMesh>> tau() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Override the pressure equation to add the mass transfer rate
            //  linearised in the pressure
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the phase and liquid's energy equations to add
            //  additional latent heat and quenching heat transfer terms
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
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

#include "phaseSurfaceBoiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
