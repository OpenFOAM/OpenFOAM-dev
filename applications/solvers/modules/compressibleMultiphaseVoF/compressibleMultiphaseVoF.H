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
    Foam::solvers::compressibleMultiphaseVoF

Description
    Solver module for the solution of multiple compressible, isothermal
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
    compressibleMultiphaseVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::multiphaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef compressibleMultiphaseVoF_H
#define compressibleMultiphaseVoF_H

#include "multiphaseVoFSolver.H"
#include "compressibleMultiphaseVoFMixture.H"
#include "compressibleMomentumTransportModel.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class compressibleMultiphaseVoF Declaration
\*---------------------------------------------------------------------------*/

class compressibleMultiphaseVoF
:
    public multiphaseVoFSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        compressibleMultiphaseVoFMixture& mixture;

        //- Reference to the phases
        UPtrListDictionary<compressibleVoFphase>& phases;


    // Thermophysical properties

        //- Reference to the mixture static pressure field
        volScalarField& p;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;

        //- Minimum pressure
        dimensionedScalar pMin;


    // Kinematic properties

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K;


    // Momentum transport

        //- Momentum transport model pointer
        autoPtr<compressible::momentumTransportModel> momentumTransport_;

        //- Momentum transport model
        compressible::momentumTransportModel& momentumTransport;


    // Cached temporary fields

        //- Continuity error
        tmp<volScalarField::Internal> contErr;


    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The flow is incompressible if all phases are incompressible
        virtual bool incompressible() const
        {
            return mixture.incompressible();
        }

        //- The flow is divergent if it is not incompressible
        //  Mass sources are not currently supported
        virtual bool divergent() const
        {
            return !incompressible();
        }

        //- Return the mixture compressibility/density
        //  Used by CorrectPhi for compressible mixtures following mesh change
        virtual tmp<volScalarField> psiByRho() const
        {
            return mixture.psiByRho();
        }

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U);

        //- Solve for the phase-fractions
        virtual void alphaPredictor();


private:

    // Private Member Functions

        //- Solve for the phase-fractions
        void alphaSolve(const dictionary& alphaControls);


public:

    //- Runtime type information
    TypeName("compressibleMultiphaseVoF");


    // Constructors

        //- Construct from region mesh
        compressibleMultiphaseVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        compressibleMultiphaseVoF
        (
            const compressibleMultiphaseVoF&
        ) = delete;


    //- Destructor
    virtual ~compressibleMultiphaseVoF();


    // Member Functions

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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleMultiphaseVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
