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
    Foam::solvers::solidDisplacement

Description
    Solver module for steady or transient segregated finite-volume solution of
    linear-elastic, small-strain deformation of a solid body, with optional
    thermal diffusion and thermal stresses.

    Solves for the displacement vector field D, also generating the stress
    tensor field sigma, including the thermal stress contribution if selected.

SourceFiles
    solidDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef solidDisplacement_H
#define solidDisplacement_H

#include "solid.H"
#include "solidDisplacementThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class solidDisplacement Declaration
\*---------------------------------------------------------------------------*/

class solidDisplacement
:
    public solid
{

protected:

    // Thermophysical properties

        solidDisplacementThermo& thermo;


    // Solution algorithm controls

        //- Switch for normal stress discretisation (required)
        //      The compact form uses snGrad (more convergent, less accurate)
        //      The non-compact form uses grad (less convergent, more accurate)
        Switch compactNormalStress;

        //- Maximum number of displacement/stress correctors per time-step
        //  Defaults to 1
        int nCorr;

        //- Convergence tolerance for the displacement/stress correctors
        //  Defaults to 0
        scalar convergenceTolerance;

        //- Acceleration factor for faster steady-state simulations
        //  Applies over-relaxation to the displacement field after solution
        //  Typically between 1 and 2, the optimal value is case dependent
        //  and should be selected to maximise convergence rate while
        //  maintaining stability.
        //  Should be set lower for higher Poisson's ratios for stability.
        scalar accFac;


    // Kinematic properties

        //- Displacement field
        volVectorField D;

        //- Young's modulus
        const volScalarField& E;

        //- Poisson's ratio
        const volScalarField& nu;

        //- Lame's coefficient
        const volScalarField mu;

        //- Lame's coefficient
        const volScalarField lambda;

        const volScalarField threeK;

        const volScalarField threeKalpha;


    // Cached temporary fields

        //- Stress field
        volSymmTensorField sigmaD;

        //- Divergence of the explicit part of the stress
        volVectorField divSigmaExp;


    // Protected Member Functions

        //- Read controls
        virtual void readControls();


public:

    //- Runtime type information
    TypeName("solidDisplacement");


    // Constructors

        //- Construct from region mesh
        solidDisplacement(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        solidDisplacement(const solidDisplacement&) = delete;


    //- Destructor
    virtual ~solidDisplacement();


    // Member Functions

        //- Called at the beginning of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the displacement equation to obtain the stress
        virtual void pressureCorrector();

        //- Correct the thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidDisplacement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
