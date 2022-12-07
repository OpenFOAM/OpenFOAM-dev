/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Solver module for for 2 compressible, non-isothermal immiscible fluids
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
    Foam::solvers::fluidSolver
    Foam::solvers::incompressibleFluid

\*---------------------------------------------------------------------------*/

#ifndef compressibleVoF_H
#define compressibleVoF_H

#include "fluidSolver.H"
#include "compressibleInterPhaseTransportModel.H"
#include "compressibleInterPhaseThermophysicalTransportModel.H"
#include "buoyancy.H"
#include "pressureReference.H"

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
    public fluidSolver
{

protected:

    // Kinematic properties

        //- Velocity field
        volVectorField U;

        //- Volumetric-flux field
        surfaceScalarField phi;


    // Phase properties

        //- The compressible two-phase mixture
        immiscibleCompressibleTwoPhaseMixture mixture;

        //- Reference to the primary phase-fraction
        volScalarField& alpha1;

        //- Switch indicating if this is a restart
        bool alphaRestart;

        scalar alphaCoNum;


    // Thermophysical properties

        //- Buoyancy force
        solvers::buoyancy buoyancy;

        //- Reference to the buoyant pressure for buoyant cases
        //  otherwise to the pressure
        volScalarField& p_rgh;

        //- Reference to the mixture continuity density field
        const volScalarField& rho;

        //- Compressibility source
        volScalarField::Internal dgdt;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference;

        //- Minimum pressure
        dimensionedScalar pMin;


    // Kinematic properties

        //- Mass flux field
        surfaceScalarField rhoPhi;

        // Phase-1 volumetric flux
        surfaceScalarField alphaPhi1;

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K;


    // Momentum transport

        //- Momentum transport model
        compressibleInterPhaseTransportModel momentumTransport;


    // Thermophysical transport

        compressibleInterPhaseThermophysicalTransportModel
            thermophysicalTransport;


    // Optional models

        //- MRF zone list
        IOMRFZoneList MRF;


    // Cached temporary fields

        tmp<volScalarField> rAU;

        tmp<volScalarField::Internal> contErr;

        //- MULES Correction
        tmp<surfaceScalarField> talphaPhi1Corr0;

        //- Pointer to the surface momentum field
        //  used to recreate the flux after mesh-change
        autoPtr<surfaceVectorField> Uf;

        //- Pointer to the momentum divergence field
        //  used in correctPhi to ensure the corrected phi has the
        //  same divergence
        autoPtr<volScalarField> divU;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;


private:

    // Private Member Functions

        //- Set rDeltaT for LTS
        virtual void setRDeltaT();

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Solve for the phase-fractions
        void alphaPredictor();


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

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual bool moveMesh();

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

        //- Correct the momentum transport modelling
        //  Newtonian, non-Newtonian or turbulent
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport modelling
        virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


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
