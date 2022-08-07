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
    Foam::solvers::isothermalFluid

Description
    Solver module for steady or transient turbulent flow of compressible
    isothermal fluids with optional mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    isothermalFluid.C

\*---------------------------------------------------------------------------*/

#ifndef isothermalFluid_H
#define isothermalFluid_H

#include "solver.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "buoyancy.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class isothermalFluid Declaration
\*---------------------------------------------------------------------------*/

class isothermalFluid
:
    public solver
{

protected:

    // Control parameters

        //- Maximum allowed Courant number
        scalar maxCo;

        //- Maximum time-step
        scalar maxDeltaT_;

        //- Switch to correct the flux after mesh change
        bool correctPhi;

        //- Switch to check the mesh Courant number after mesh change
        bool checkMeshCourantNo;


    // Thermophysical properties

        //- Pointer to the fluid thermophysical properties
        autoPtr<fluidThermo> thermo_;

        //- Reference to the fluid thermophysical properties
        fluidThermo& thermo;

        //- Reference to the pressure field
        volScalarField& p;

        //- The continuity density field
        volScalarField rho;

        //- Rate of change of the pressure
        //  Used in the enthalpy equation
        volScalarField dpdt;


    // Optional buoyancy

        //- Pointer to the optional buoyancy force
        //  Case is considered buoyant if the p_rgh field exists
        autoPtr<solvers::buoyancy> buoyancy;

        //- Reference to the buoyant pressure for buoyant cases
        //  otherwise to the pressure
        volScalarField& p_rgh;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference;


    // Kinematic properties

        //- Velocity field
        volVectorField U;

        //- Mass-flux field
        surfaceScalarField phi;

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K;


    // Momentum transport

        //- Pointer to the momentum transport model
        autoPtr<compressible::momentumTransportModel> turbulence;


    // Continuity properties

        //- Initial mass in the region
        dimensionedScalar initialMass;

        //- Current cumulative continuity error
        scalar cumulativeContErr;


    // Optional models

        //- MRF zone list
        IOMRFZoneList MRF;


    // Time-step controls

        //- Current maximum Courant number
        scalar CoNum;

        //- Current mean Courant number
        scalar meanCoNum;


    // Cached temporary fields

        //- Pointer to the surface momentum field
        //  used to recreate the flux after mesh-change
        autoPtr<surfaceVectorField> rhoUf;

        //- Pointer to the vol momentum field
        //  used for mesh-change to set rhoUf for introduced faces
        autoPtr<volVectorField> rhoU;

        //- Pointer to the vol momentum divergence field
        //  used in correctPhi to ensure the corrected phi has the
        //  same divergence
        autoPtr<volScalarField> divrhoU;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;


private:

    // Private Member Data

        //- Pointer to the demand driven fvModels MeshObject
        mutable Foam::fvModels* fvModelsPtr;

        //- Pointer to the demand driven fvConstraints MeshObject
        mutable Foam::fvConstraints* fvConstraintsPtr;


    // Private Member Functions

        //- Read controls
        void read();

        //- Set rDeltaT for LTS
        virtual void setRDeltaT();

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Check mesh Courant numbers for moving mesh cases
        void meshCourantNo() const;

        //- Calculate and print the continuity errors
        void continuityErrors();

        //- Construct the continuity equation and correct the density
        void correctDensity();

        //- Construct the pressure equation
        //  and correct the pressure and velocity
        void correctPressure();

        //- Construct the buoyant pressure equation
        //  and correct the pressure and velocity
        void correctBuoyantPressure();


public:

    //- Runtime type information
    TypeName("isothermalFluid");


    // Constructors

        //- Construct from region mesh and thermophysical properties
        isothermalFluid(fvMesh& mesh, autoPtr<fluidThermo>);

        //- Construct from region mesh
        isothermalFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        isothermalFluid(const isothermalFluid&) = delete;


    //- Destructor
    virtual ~isothermalFluid();


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
        void operator=(const isothermalFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
