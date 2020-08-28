/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::SmoothSolver

Description
    Iterative solver for symmetric and asymmetric matrices which uses a
    run-time selected smoother e.g. GaussSeidel to converge the solution to
    the required tolerance.  To improve efficiency, the residual is evaluated
    after every nSweeps smoothing iterations.

SourceFiles
    SmoothSolver.C

\*---------------------------------------------------------------------------*/

#ifndef SmoothSolver_H
#define SmoothSolver_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SmoothSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class SmoothSolver
:
    public LduMatrix<Type, DType, LUType>::solver
{

protected:

    // Protected data

        //- Number of sweeps before the evaluation of residual
        label nSweeps_;

        //- Read the control parameters from the controlDict_
        virtual void readControls();


public:

    //- Runtime type information
    TypeName("SmoothSolver");


    // Constructors

        //- Construct from matrix components and solver data dictionary
        SmoothSolver
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverDict
        );


    // Member Functions

        //- Solve the matrix with this solver
        virtual SolverPerformance<Type> solve(Field<Type>& psi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SmoothSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
