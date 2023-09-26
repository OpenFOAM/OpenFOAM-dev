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
    Foam::solvers::incompressibleDriftFlux

Description
    Solver module for 2 incompressible fluids using the mixture approach with
    the drift-flux approximation for relative motion of the phases, with
    optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved with mixture transport modelling in which a
    single laminar, RAS or LES model is selected to model the momentum stress.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    incompressibleDriftFlux.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseSolver
    Foam::solvers::compressibleVoF

\*---------------------------------------------------------------------------*/

#ifndef incompressibleDriftFlux_H
#define incompressibleDriftFlux_H

#include "twoPhaseSolver.H"
#include "incompressibleDriftFluxMixture.H"
#include "relativeVelocityModel.H"
#include "momentumTransportModel.H"
#include "compressibleMomentumTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleDriftFlux Declaration
\*---------------------------------------------------------------------------*/

class incompressibleDriftFlux
:
    public twoPhaseSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        incompressibleDriftFluxMixture& mixture;


    // Thermophysical properties

        //- Static pressure field
        volScalarField p;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Momentum transport

        //- Pointer to the dispersed phase relative velocity model
        autoPtr<relativeVelocityModel> relativeVelocity;

        //- Pointer to the momentum transport model
        autoPtr<compressible::momentumTransportModel> momentumTransport;


private:

    // Private Member Functions

        //- Correct the cached Courant numbers
        virtual void correctCoNum();


protected:

    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The flow is incompressible
        virtual bool incompressible() const
        {
            return true;
        }

        //- Is the flow divergent?
        //  i.e. includes phase-fraction sources
        virtual bool divergent() const
        {
            return
                fvModels().addsSupToField(alpha1.name())
             || fvModels().addsSupToField(alpha2.name());
        }

        //- Return the mixture compressibility/density
        //  Not required for incompressible fluids
        virtual tmp<volScalarField> psiByRho() const
        {
            return tmp<volScalarField>(nullptr);
        }

        virtual tmp<surfaceScalarField> alphaPhi
        (
            const surfaceScalarField& phi,
            const volScalarField& alpha,
            const dictionary& alphaControls
        );

        //- Adjust the rDeltaT in the vicinity of the interface
        virtual void setInterfaceRDeltaT(volScalarField& rDeltaT);

        //- Calculate the alpha equation sources
        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp
        );

        //- Correct the interface properties following mesh-change
        //  and phase-fraction update
        virtual void correctInterface();

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U)
        {
            return
                relativeVelocity->divDevTau()
              + momentumTransport->divDevTau(U);
        }


public:

    //- Runtime type information
    TypeName("incompressibleDriftFlux");


    // Constructors

        //- Construct from region mesh
        incompressibleDriftFlux(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleDriftFlux(const incompressibleDriftFlux&) = delete;


    //- Destructor
    virtual ~incompressibleDriftFlux();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

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
        void operator=(const incompressibleDriftFlux&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
