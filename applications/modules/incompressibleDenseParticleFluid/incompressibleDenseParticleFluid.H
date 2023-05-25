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
    Foam::solvers::incompressibleDenseParticleFluid

Description

    Solver module for transient flow of incompressible isothermal fluids coupled
    with particle clouds including the effect of the volume fraction of
    particles on the continuous phase, with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, constraining or limiting
    the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    incompressibleDenseParticleFluid.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::incompressibleFluid

\*---------------------------------------------------------------------------*/

#ifndef incompressibleDenseParticleFluid_H
#define incompressibleDenseParticleFluid_H

#include "fluidSolver.H"
#include "viscosityModel.H"
#include "phaseIncompressibleMomentumTransportModel.H"
#include "pressureReference.H"
#include "uniformDimensionedFields.H"
#include "parcelClouds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleDenseParticleFluid Declaration
\*---------------------------------------------------------------------------*/

class incompressibleDenseParticleFluid
:
    public fluidSolver
{

protected:

    //- Name of the continuous phase
    const word continuousPhaseName;


    // Pressure

        //- Pressure field
        volScalarField p_;

        //- Pressure reference
        Foam::pressureReference pressureReference;


    // Buoyancy

        //- Acceleration due to gravity
        uniformDimensionedVectorField g;


    // Kinematic properties

        //- Continuous phase velocity field
        volVectorField Uc_;

        //- Continuous phase flux field
        surfaceScalarField phic_;

        //- Continuous phase kinematic viscosity model
        autoPtr<viscosityModel> viscosity;


    // Continuous phase properties for Lagrangian cloud

        //- Continuous phase density
        volScalarField rhoc;

        //- Continuous phase viscosity
        volScalarField muc;

        //- Continuous phase-fraction
        volScalarField alphac_;

        //- Minimum continuous phase-fraction
        scalar alphacMin;

        //- Interpolated continuous phase-fraction
        surfaceScalarField alphacf;

        //- Continuous phase volumetric-flux field
        surfaceScalarField alphaPhic;


    // Momentum transport

        //- Pointer to the momentum transport model
        autoPtr<phaseIncompressible::momentumTransportModel> momentumTransport;


    // Dispersed phase Lagrangian clouds

        parcelClouds& clouds;


    // Cached temporary fields

        //- Dispersed phase drag force
        autoPtr<volVectorField> Fd;

        // Continuous-dispersed phase drag coefficient
        autoPtr<volScalarField> Dc;

        //- Face continuous-dispersed phase drag coefficient
        autoPtr<surfaceScalarField> Dcf;

        //- Effective volumetric flux of the dispersed phase
        autoPtr<surfaceScalarField> phid;

        //- Pointer to the surface momentum field
        //  used to recreate the flux after mesh-change
        autoPtr<surfaceVectorField> Ucf;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUcEqn;


    // Protected Member Functions

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Calculate and print the continuity errors
        void continuityErrors();

        //- Construct the pressure equation
        //  and correct the pressure and velocity
        virtual void correctPressure();


public:

    // Public Data

        //- Reference to the pressure field
        const volScalarField& p;

        //- Reference to the continuous phase velocity field
        const volVectorField& Uc;

        //- Reference to the continuous phase volumetric-flux field
        const surfaceScalarField& phic;

        //- Reference continuous phase-fraction
        const volScalarField& alphac;


    //- Runtime type information
    TypeName("incompressibleDenseParticleFluid");


    // Constructors

        //- Construct from region mesh
        incompressibleDenseParticleFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleDenseParticleFluid
        (
            const incompressibleDenseParticleFluid&
        ) = delete;


    //- Destructor
    virtual ~incompressibleDenseParticleFluid();


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
        void operator=(const incompressibleDenseParticleFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
