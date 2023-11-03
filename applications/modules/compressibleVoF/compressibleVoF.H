/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::solvers::compressibleVoF

Description
    Solver module for 2 compressible, non-isothermal immiscible fluids
    using a VOF (volume of fluid) phase-fraction based interface capturing
    approach, with optional mesh motion and mesh topology changes including
    adaptive re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Either mixture or two-phase transport modelling may be selected.  In the
    mixture approach a single laminar, RAS or LES model is selected to model the
    momentum stress.  In the Euler-Euler two-phase approach separate laminar,
    RAS or LES selected models are selected for each of the phases.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    compressibleVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseVoFSolver
    Foam::solvers::incompressibleVoF

\*---------------------------------------------------------------------------*/

#ifndef compressibleVoF_H
#define compressibleVoF_H

#include "twoPhaseVoFSolver.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "compressibleInterPhaseTransportModel.H"
#include "compressibleInterPhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class compressibleVoF Declaration
\*---------------------------------------------------------------------------*/

class compressibleVoF
:
    public twoPhaseVoFSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        compressibleTwoPhaseVoFMixture& mixture_;


    // Thermophysical properties

        //- Reference to the mixture static pressure field
        volScalarField& p;

        //- Compressibility source
        volScalarField::Internal vDot;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Kinematic properties

        // Phase-1 mass-flux
        surfaceScalarField alphaRhoPhi1;

        // Phase-2 mass-flux
        surfaceScalarField alphaRhoPhi2;

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K;


    // Momentum transport

        //- Momentum transport model
        compressibleInterPhaseTransportModel momentumTransport;


    // Thermophysical transport

        //- Thermophysical transport model
        compressibleInterPhaseThermophysicalTransportModel
            thermophysicalTransport;


    // Cached temporary fields

        //- Phase-1 continuity error
        tmp<volScalarField::Internal> contErr1;

        //- Phase-2 continuity error
        tmp<volScalarField::Internal> contErr2;


    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The fluid is incompressible if both phases are incompressible
        virtual bool incompressible() const
        {
            return mixture.incompressible();
        }

        //- Compressible flow is divergent
        virtual bool divergent() const
        {
            return
               !incompressible()
             || fvModels().addsSupToField(alpha1.name())
             || fvModels().addsSupToField(alpha2.name());
        }

        //- Return the mixture compressibility/density
        //  Used by CorrectPhi for compressible mixtures following mesh change
        virtual tmp<volScalarField> psiByRho() const
        {
            return mixture.psiByRho();
        }

        //- Calculate the alpha equation sources
        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp,
            const dictionary& alphaControls
        );

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U);

public:

    // Access to the physical state

        //- The compressible two-phase mixture
        const compressibleTwoPhaseVoFMixture& mixture;


public:

    //- Runtime type information
    TypeName("compressibleVoF");


    // Constructors

        //- Construct from region mesh
        compressibleVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        compressibleVoF(const compressibleVoF&) = delete;


    //- Destructor
    virtual ~compressibleVoF();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
