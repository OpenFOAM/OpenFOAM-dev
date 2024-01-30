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
    Foam::solvers::VoFSolver

Description
    Base solver module base-class for the solution of immiscible fluids using
    a VOF (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    VoFSolver.C

See also
    Foam::solvers::fluidSolver

\*---------------------------------------------------------------------------*/

#ifndef VoFSolver_H
#define VoFSolver_H

#include "fluidSolver.H"
#include "VoFMixture.H"
#include "buoyancy.H"
#include "pressureReference.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class VoFSolver Declaration
\*---------------------------------------------------------------------------*/

class VoFSolver
:
    public fluidSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        autoPtr<VoFMixture> mixturePtr_;

        VoFMixture& mixture_;

        //- Name of the alpha convection scheme
        const word divAlphaName;

        //- Phase-fraction flux Courant number
        scalar alphaCoNum;


    // Kinematic properties

        //- Velocity field
        volVectorField U_;

        //- Volumetric flux field
        surfaceScalarField phi_;


    // Thermophysical properties

        //- Buoyancy force
        solvers::buoyancy buoyancy;

        //- Reference to the buoyant pressure for buoyant cases
        //  otherwise to the pressure
        volScalarField& p_rgh;

        //- Reference to the mixture continuity density field
        const volScalarField& rho;


    // Kinematic properties

        //- Mass flux field
        surfaceScalarField rhoPhi;


    // Optional models

        //- MRF zone list
        IOMRFZoneList MRF;


    // Cached temporary fields

        //- Inverse momentum equation diagonal
        //  Used to correct phi following mesh changes
        autoPtr<volScalarField> rAU;

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
        void setRDeltaT();


protected:

    // Protected Member Functions

        //- Set or update the cached rAU
        void setrAU(const fvVectorMatrix& UEqn);

        //- Clear the cached rAU is no longer needed
        void clearrAU();

        //- Correct the cached Courant numbers
        virtual void correctCoNum() = 0;

        //- Adjust the rDeltaT in the vicinity of the interface
        virtual void setInterfaceRDeltaT(volScalarField& rDeltaT) = 0;

        //- Calculate and print the continuity errors
        void continuityErrors();

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const = 0;

        //- Is the flow incompressible?
        virtual bool incompressible() const = 0;

        //- Is the flow divergent?
        //  i.e. compressible or include phase-fraction sources
        virtual bool divergent() const = 0;

        //- Return the mixture compressibility/density
        //  Used by CorrectPhi for compressible mixtures following mesh change
        virtual tmp<volScalarField> psiByRho() const = 0;

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) = 0;

        //- Correct the interface properties following mesh-change
        //  and phase-fraction update
        virtual void correctInterface() = 0;

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const = 0;


public:

    // Public Data

        //- Reference to the VoF mixture
        const VoFMixture& mixture;

        //- Reference to the velocity field
        const volVectorField& U;

        //- Reference to the mass-flux field
        const surfaceScalarField& phi;


    //- Runtime type information
    TypeName("VoFSolver");


    // Constructors

        //- Construct from region mesh
        VoFSolver(fvMesh& mesh, autoPtr<VoFMixture>);

        //- Disallow default bitwise copy construction
        VoFSolver(const VoFSolver&) = delete;


    //- Destructor
    virtual ~VoFSolver();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor() = 0;

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor() = 0;

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector() = 0;

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
