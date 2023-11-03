/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::solvers::twoPhaseSolver

Description
    Solver module base-class for 2 immiscible fluids, with optional
    mesh motion and mesh topology changes including adaptive re-meshing.

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
    twoPhaseSolver.C

See also
    Foam::solvers::fluidSolver

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseSolver_H
#define twoPhaseSolver_H

#include "VoFSolver.H"
#include "twoPhaseVoFMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class twoPhaseSolver Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSolver
:
    public VoFSolver
{

protected:

    // Phase properties

        //- Reference to the twoPhaseVoFMixture
        twoPhaseVoFMixture& mixture;

        //- Reference to the phase1-fraction
        volScalarField& alpha1;

        //- Reference to the phase2-fraction
        volScalarField& alpha2;

        //- Switch indicating if this is a restart
        bool alphaRestart;


    // Kinematic properties

        // Phase-1 volumetric flux
        surfaceScalarField alphaPhi1;


    // Cached temporary fields

        //- MULES Correction
        tmp<surfaceScalarField> talphaPhi1Corr0;


private:

    // Private Member Functions

        //- Solve for the phase-fractions
        void alphaSolve(const dictionary& alphaControls);


protected:

    // Protected Member Functions

        virtual tmp<surfaceScalarField> alphaPhi
        (
            const surfaceScalarField& phi,
            const volScalarField& alpha,
            const dictionary& alphaControls
        );

        //- Solve for the phase-fractions
        void alphaPredictor();

        //- Calculate the alpha equation sources
        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp,
            const dictionary& alphaControls
        ) = 0;

        //- Correct the interface properties following mesh-change
        //  and phase-fraction update
        virtual void correctInterface() = 0;

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const = 0;

        //- Construct and solve the incompressible pressure equation
        //  in the PISO loop
        void incompressiblePressureCorrector(volScalarField& p);


public:

    //- Runtime type information
    TypeName("twoPhaseSolver");


    // Constructors

        //- Construct from region mesh
        twoPhaseSolver(fvMesh& mesh, autoPtr<twoPhaseVoFMixture>);

        //- Disallow default bitwise copy construction
        twoPhaseSolver(const twoPhaseSolver&) = delete;


    //- Destructor
    virtual ~twoPhaseSolver();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const twoPhaseSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
