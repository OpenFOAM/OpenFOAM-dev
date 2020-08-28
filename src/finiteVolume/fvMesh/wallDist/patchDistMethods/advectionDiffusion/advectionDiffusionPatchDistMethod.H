/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::patchDistMethods::advectionDiffusion

Description
    Calculation of approximate distance to nearest patch for all cells and
    boundary by solving the Eikonal equation in advection form with diffusion
    smoothing.

    If the diffusion coefficient is set to 0 this method is exact in principle
    but the numerical schemes used are not rendering the scheme approximate, but
    more accurate than the Poisson method.  Also many models relying on the
    distance to the wall benefit from this field being smooth and monotonic so
    the addition of diffusion smoothing improves both the convergence and
    stability of the solution of the Eikonal equation and the behavior of the
    models using the distance field generated.  However, it is not clear what
    the optimum value for the diffusion coefficient epsilon should be; a
    default value of 0.1 is provided but higher values may be preferable under
    some circumstances.

    Convergence is accelerated by first generating an approximate solution
    using one of the simpler methods, e.g. Poisson.  The method used for
    this prediction step is specified in the 'advectionDiffusionCoeffs'
    sub-dictionary, see below.

    References:
    \verbatim
        P.G. Tucker, C.L. Rumsey, R.E. Bartels, R.T. Biedron,
        "Transport equation based wall distance computations aimed at flows
         with time-dependent geometry",
        NASA/TM-2003-212680, December 2003.
    \endverbatim

    Example of the wallDist specification in fvSchemes:
    \verbatim
        laplacianSchemes
        {
            .
            .
            laplacian(yPsi) Gauss linear corrected;
            laplacian(yWall) Gauss linear corrected;
            .
            .
        }

        wallDist
        {
            method advectionDiffusion;

            // Optional entry enabling the calculation
            // of the normal-to-wall field
            nRequired false;

            advectionDiffusionCoeffs
            {
                method    Poisson;
                epsilon   0.1;
                tolerance 1e-3;
                maxIter   10;
            }
        }
    \endverbatim
    Also the solver specification for yWall is required in fvSolution, e.g.
    for simple cases:
    \verbatim
        yPsi
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-4;
            relTol          0;
        }

        yWall
        {
            solver          PBiCG;
            preconditioner  DILU;
            tolerance       1e-4;
            relTol          0;
        }

        relaxationFactors
        {
            equations
            {
                .
                .
                yWall           1;
                .
                .
            }
        }

    or for more complex cases:

        yPsi
        {
            solver          GAMG;
            smoother        GaussSeidel;
            cacheAgglomeration true;
            nCellsInCoarsestLevel 10;
            agglomerator    faceAreaPair;
            mergeLevels     1;
            tolerance       1e-4;
            relTol          0;
        }

        yWall
        {
            solver          GAMG;
            smoother        symGaussSeidel;
            cacheAgglomeration true;
            nCellsInCoarsestLevel 10;
            agglomerator    faceAreaPair;
            mergeLevels     1;
            tolerance       1e-4;
            relTol          0;
        }

        relaxationFactors
        {
            equations
            {
                .
                .
                yWall           0.7;
                .
                .
            }
        }
    \endverbatim

See also
    Foam::patchDistMethod::Poisson
    Foam::patchDistMethod::meshWave
    Foam::wallDist

SourceFiles
    advectionDiffusionPatchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef advectionDiffusionPatchDistMethod_H
#define advectionDiffusionPatchDistMethod_H

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{

/*---------------------------------------------------------------------------*\
                     Class advectionDiffusion Declaration
\*---------------------------------------------------------------------------*/

class advectionDiffusion
:
    public patchDistMethod
{
    // Private Member Data

        //- Sub-dictionary of coefficients
        dictionary coeffs_;

        //- Run-time selected method to predict the distance-to-wall field
        autoPtr<patchDistMethod> pdmPredictor_;

        //- Diffusion coefficient multiplying y*laplacian(y)
        scalar epsilon_;

        //- Convergence tolerance for the iterations of the advection-diffusion
        //  equation to correct the distance-to-patch and normal-to-patch fields
        scalar tolerance_;

        //- Maximum number of iterations of the advection-diffusion equation
        //  to correct the distance-to-patch and normal-to-patch fields
        int maxIter_;

        //- Flag to indicate the predictor step has been executed
        bool predicted_;


public:

    //- Runtime type information
    TypeName("advectionDiffusion");


    // Constructors

        //- Construct from coefficients dictionary, mesh
        //  and fixed-value patch set
        advectionDiffusion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Disallow default bitwise copy construction
        advectionDiffusion(const advectionDiffusion&) = delete;


    // Member Functions

        //- Correct the given distance-to-patch field
        virtual bool correct(volScalarField& y);

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(volScalarField& y, volVectorField& n);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const advectionDiffusion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
