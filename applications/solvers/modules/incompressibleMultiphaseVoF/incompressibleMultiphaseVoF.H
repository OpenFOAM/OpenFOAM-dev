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
    Foam::solvers::incompressibleMultiphaseVoF

Description
    Solver module for the solution of multiple incompressible, isothermal
    immiscible fluids using a VOF (volume of fluid) phase-fraction based
    interface capturing approach, with optional mesh motion and mesh topology
    changes including adaptive re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    A mixture approach for momentum transport is provided in which a single
    laminar, RAS or LES model is selected to model the momentum stress.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

SourceFiles
    incompressibleMultiphaseVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::multiphaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef incompressibleMultiphaseVoF_H
#define incompressibleMultiphaseVoF_H

#include "multiphaseVoFSolver.H"
#include "incompressibleMultiphaseVoFMixture.H"
#include "incompressibleMomentumTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleMultiphaseVoF Declaration
\*---------------------------------------------------------------------------*/

class incompressibleMultiphaseVoF
:
    public multiphaseVoFSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        incompressibleMultiphaseVoFMixture& mixture;

        //- Reference to the phases
        UPtrListDictionary<incompressibleVoFphase>& phases;


    // Thermophysical properties

        //- Static pressure field
        volScalarField p;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Momentum transport

        //- Momentum transport model pointer
        autoPtr<incompressible::momentumTransportModel> momentumTransport_;

        //- Momentum transport model
        incompressible::momentumTransportModel& momentumTransport;


    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- Is the flow divergent?
        //  i.e. includes phase-fraction sources
        virtual bool divergent()
        {
            return false;
        }

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U)
        {
            return momentumTransport.divDevTau(rho, U);
        }

        //- Solve for the phase-fractions
        virtual void alphaPredictor();


private:

    // Private Member Functions

        //- Solve for the phase-fractions
        void alphaSolve(const dictionary& alphaControls);


public:

    //- Runtime type information
    TypeName("incompressibleMultiphaseVoF");


    // Constructors

        //- Construct from region mesh
        incompressibleMultiphaseVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleMultiphaseVoF
        (
            const incompressibleMultiphaseVoF&
        ) = delete;


    //- Destructor
    virtual ~incompressibleMultiphaseVoF();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

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
        void operator=(const incompressibleMultiphaseVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
