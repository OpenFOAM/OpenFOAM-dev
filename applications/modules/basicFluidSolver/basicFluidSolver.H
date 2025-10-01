/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::solvers::basicFluidSolver

Description
    Base solver module for fluid solvers.

    Provides Courant number time-step control and continuity checking.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    basicFluidSolver.C

\*---------------------------------------------------------------------------*/

#ifndef basicFluidSolver_H
#define basicFluidSolver_H

#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class basicFluidSolver Declaration
\*---------------------------------------------------------------------------*/

class basicFluidSolver
:
    public solver
{
    // Control parameters

        //- Maximum allowed Courant number
        scalar maxCo;

        //- Maximum time-step
        scalar maxDeltaT_;

        //- Switch to check the mesh Courant number after mesh change
        bool checkMeshCourantNo;


    // Continuity properties

        //- Current cumulative continuity error
        scalar cumulativeContErr;


    // Private Member Functions

        //- Correct the cached Courant numbers
        template<class RhoType>
        inline void correctCoNum
        (
            const RhoType& rho,
            const surfaceScalarField& phi
        );


protected:

        //- Return true if the solver's dependencies have been modified
        virtual bool dependenciesModified() const;

        //- Read controls
        virtual bool read();

        //- Switch to correct the flux after mesh change
        bool correctPhi;

        //- Current maximum Courant number for time-step control
        scalar CoNum_;

        //- Check mesh Courant numbers for moving mesh cases
        void meshCourantNo() const;

        //- Correct the cached Courant numbers
        void correctCoNum(const surfaceScalarField& phi);

        //- Correct the cached Courant numbers
        void correctCoNum
        (
            const volScalarField& rho,
            const surfaceScalarField& phi
        );

        //- Calculate and print the continuity errors
        void continuityErrors
        (
            const surfaceScalarField& phi
        );

        //- Calculate and print the continuity errors
        void continuityErrors
        (
            const volScalarField& rho,
            const volScalarField& thermoRho,
            const surfaceScalarField& phi
        );


public:

    // Public Data

        //- Current maximum Courant number for time-step control
        const scalar& CoNum;


public:

    //- Runtime type information
    TypeName("basicFluidSolver");


    // Constructors

        //- Construct from region mesh
        basicFluidSolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        basicFluidSolver(const basicFluidSolver&) = delete;


    //- Destructor
    virtual ~basicFluidSolver();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicFluidSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
