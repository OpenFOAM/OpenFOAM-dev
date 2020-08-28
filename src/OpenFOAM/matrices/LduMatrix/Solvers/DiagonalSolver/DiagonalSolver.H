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
    Foam::DiagonalSolver

Description
    Foam::DiagonalSolver

SourceFiles
    DiagonalSolver.C

\*---------------------------------------------------------------------------*/

#ifndef DiagonalSolver_H
#define DiagonalSolver_H

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DiagonalSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class DiagonalSolver
:
    public LduMatrix<Type, DType, LUType>::solver
{

public:

    //- Runtime type information
    TypeName("diagonal");


    // Constructors

        //- Construct from matrix
        DiagonalSolver
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverDict
        );

        //- Disallow default bitwise copy construction
        DiagonalSolver(const DiagonalSolver&) = delete;


    // Member Functions

        //- Read and reset the solver parameters from the given dictionary
        void read(const dictionary& solverDict);

        //- Solve the matrix with this solver
        virtual SolverPerformance<Type> solve(Field<Type>& psi) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DiagonalSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagonalSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
