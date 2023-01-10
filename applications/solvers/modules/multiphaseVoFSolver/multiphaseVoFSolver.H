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
    Foam::solvers::multiphaseVoFSolver

Description
    Base solver module for the solution of multiple immiscible fluids using a
    VOF (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

SourceFiles
    multiphaseVoFSolver.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::incompressibleMultiphaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef multiphaseVoFSolver_H
#define multiphaseVoFSolver_H

#include "VoFSolver.H"
#include "multiphaseVoFMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseVoFSolver Declaration
\*---------------------------------------------------------------------------*/

class multiphaseVoFSolver
:
    public VoFSolver
{

protected:

    // Phase properties

        //- Reference to the multiphaseVoFMixture
        multiphaseVoFMixture& mixture;

        //- Reference to the phases
        PtrListDictionary<VoFphase>& phases;


private:

    // Private Member Functions

        //- Correct the cached Courant numbers
        void correctCoNum();


protected:

    // Protected Member Functions

        //- Adjust the rDeltaT in the vicinity of the interface
        virtual void setInterfaceRDeltaT(volScalarField& rDeltaT);

        //- Solve for the phase-fractions
        virtual void alphaPredictor() = 0;

        //- Correct the interface properties following mesh-change
        //  and phase-fraction update
        virtual void correctInterface();

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;


public:

    //- Runtime type information
    TypeName("multiphaseVoFSolver");


    // Constructors

        //- Construct from region mesh
        multiphaseVoFSolver(fvMesh& mesh, autoPtr<multiphaseVoFMixture>);

        //- Disallow default bitwise copy construction
        multiphaseVoFSolver(const multiphaseVoFSolver&) = delete;


    //- Destructor
    virtual ~multiphaseVoFSolver();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseVoFSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
