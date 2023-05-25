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
    Foam::solvers::twoPhaseVoFSolver

Description
    Solver module base-class for for 2 immiscible fluids using a VOF (volume
    of fluid) phase-fraction based interface capturing approach, with optional
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
    twoPhaseVoFSolver.C

See also
    Foam::solvers::fluidSolver

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseVoFSolver_H
#define twoPhaseVoFSolver_H

#include "twoPhaseSolver.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class twoPhaseVoFSolver Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseVoFSolver
:
    public twoPhaseSolver
{

protected:

    // Interface properties

        interfaceProperties interface;


private:

    // Private Member Functions

        //- Correct the cached Courant numbers
        virtual void correctCoNum();


protected:

    // Protected Member Functions

        //- Adjust the rDeltaT in the vicinity of the interface
        virtual void setInterfaceRDeltaT(volScalarField& rDeltaT);

        //- Correct the interface properties following mesh-change
        //  and phase-fraction update
        virtual void correctInterface();

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;


public:

    //- Runtime type information
    TypeName("twoPhaseVoFSolver");


    // Constructors

        //- Construct from region mesh
        twoPhaseVoFSolver(fvMesh& mesh, autoPtr<twoPhaseVoFMixture>);

        //- Disallow default bitwise copy construction
        twoPhaseVoFSolver(const twoPhaseVoFSolver&) = delete;


    //- Destructor
    virtual ~twoPhaseVoFSolver();


    // Member Functions


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const twoPhaseVoFSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
