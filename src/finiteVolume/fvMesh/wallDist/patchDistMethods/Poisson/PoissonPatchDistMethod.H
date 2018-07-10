/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::patchDistMethods::Poisson

Description
    Calculation of approximate distance to nearest patch for all cells and
    boundary by solving Poisson's equation.

    References:
    \verbatim
        D.B. Spalding,
        "Calculation of turbulent heat transfer in cluttered spaces",
        Proc. 10th Int. Heat Transfer Conference, Brighton, UK, (1994).

        E. Fares and W. Schroder,
        "Differential Equation for Approximate Wall Distance",
        Int.J.Numer.Meth., 39:743-762, (2002).

        P.G. Tucker,
        "Differential equation based wall distance computation for DES and
         RANS",
        J.Comp.Phys., Vol. 190, Issue 1, 1 st September, pp. 229-248 (2003)
    \endverbatim

    Example of the wallDist specification in fvSchemes:
    \verbatim
        laplacianSchemes
        {
            .
            .
            laplacian(yPsi) Gauss linear corrected;
            .
            .
        }

        wallDist
        {
            method Poisson;

            // Optional entry enabling the calculation
            // of the normal-to-wall field
            nRequired false;
        }
    \endverbatim
    Also the solver specification for yPsi is required in fvSolution, e.g.
    for simple cases:
    \verbatim
        yPsi
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-5;
            relTol          0;
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
            tolerance       1e-5;
            relTol          0;
        }
    \endverbatim

See also
    Foam::patchDistMethod::meshWave
    Foam::wallDist

SourceFiles
    PoissonPatchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef PoissonPatchDistMethod_H
#define PoissonPatchDistMethod_H

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{

/*---------------------------------------------------------------------------*\
                          Class Poisson Declaration
\*---------------------------------------------------------------------------*/

class Poisson
:
    public patchDistMethod
{
    // Private Member Data

        //- Cache yPsi for moving meshes
        tmp<volScalarField> tyPsi_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Poisson(const Poisson&);

        //- Disallow default bitwise assignment
        void operator=(const Poisson&);


public:

    //- Runtime type information
    TypeName("Poisson");


    // Constructors

        //- Construct from coefficients dictionary, mesh
        //  and fixed-value patch set
        Poisson
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Construct from mesh and fixed-value patch set
        Poisson
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );


    // Member Functions

        //- Correct the given distance-to-patch field
        virtual bool correct(volScalarField& y);

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(volScalarField& y, volVectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
