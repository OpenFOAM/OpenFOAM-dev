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
    Foam::solvers::multiphaseEuler

Description
    Solver module for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentum, heat and mass transfer.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    multiphaseEuler.C

See also
    Foam::solvers::compressibleVoF
    Foam::solvers::fluidSolver
    Foam::solvers::incompressibleFluid

\*---------------------------------------------------------------------------*/

#ifndef multiphaseEuler_H
#define multiphaseEuler_H

#include "fluidSolver.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "buoyancy.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseEuler Declaration
\*---------------------------------------------------------------------------*/

class multiphaseEuler
:
    public fluidSolver
{

protected:

    // Controls

        //- Momentum equation predictor switch
        //  Defaults to false
        Switch predictMomentum;

        //- Cell/face momentum equation switch
        //  Defaults to false, i.e. uses the cell momentum equation
        Switch faceMomentum;

        //- Cell/face drag correction for cell momentum corrector
        //  Defaults to false
        Switch dragCorrection;

        //- Number of energy correctors
        //  Used to improve stability of phase-change sibulations
        //  Defaults to 1
        int nEnergyCorrectors;


    //- Optional LTS reciprocal time-step field
    tmp<volScalarField> trDeltaT;

    //- Optional LTS reciprocal face time-step field
    tmp<surfaceScalarField> trDeltaTf;

    //- Buoyancy force
    solvers::buoyancy buoyancy;


    // Phase properties

        autoPtr<phaseSystem> fluidPtr_;

        phaseSystem& fluid_;

        phaseSystem::phaseModelList& phases_;

        phaseSystem::phaseModelPartialList& movingPhases_;

        surfaceScalarField& phi_;


    // Thermophysical properties

        volScalarField& p_;

        //- Reference to the buoyant pressure for buoyant cases
        //  otherwise to the pressure
        volScalarField& p_rgh;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference;


    // Optional models

        const IOMRFZoneList& MRF;


    // Cached temporary fields

        //- Temporary phase momentum matrices
        //  shared between the momentum predictor and pressure corrector
        PtrList<fvVectorMatrix> UEqns;

        //- Temporary storage for the reciprocal momentum equation diagonal
        //  Used by the phase-fraction predictor and pressure corrector
        PtrList<volScalarField> rAs;

        //- Stored divU from the previous mesh so that it can be
        //  mapped and used in correctPhi to ensure the corrected phi
        //  has the same divergence
        autoPtr<volScalarField> divU;

        //- Read controls
        virtual bool read();


private:

    // Private Member Functions

        //- Set rDeltaT for LTS
        virtual void setRDeltaT();

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Composition predictor
        //  called in prePredictor() after phase-fraction prediction
        void compositionPredictor();

        //- Energy and temperature predictor
        void energyPredictor();

        //- Cell-based momentum equation predictor
        void cellMomentumPredictor();

        //- Face-based momentum equation predictor
        void faceMomentumPredictor();

        //- Cell-based pressure equation predictor
        void cellPressureCorrector();

        //- Face-based pressure equation predictor
        void facePressureCorrector();

        //- Return the list of pressure equation compressibility contributions
        PtrList<fvScalarMatrix> compressibilityEqns
        (
            const PtrList<volScalarField>& dmdts,
            const PtrList<volScalarField>& d2mdtdps
        ) const;


public:

    // Public Data

        //- Reference to the multiphase fluid
        const phaseSystem& fluid;

        //- Reference to the phases
        const phaseSystem::phaseModelList& phases;

        //- Reference to the moving phases
        const phaseSystem::phaseModelPartialList& movingPhases;

        //- Reference to the pressure field
        const volScalarField& p;

        //- Reference to the mass-flux field
        const surfaceScalarField& phi;


    //- Runtime type information
    TypeName("multiphaseEuler");


    // Constructors

        //- Construct from region mesh
        multiphaseEuler(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        multiphaseEuler(const multiphaseEuler&) = delete;


    //- Destructor
    virtual ~multiphaseEuler();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

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

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseEuler&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
