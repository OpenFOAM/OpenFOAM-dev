/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::LLTMatrix

Description
    Templated class to perform the Cholesky decomposition on a
    symmetric positive-definite matrix.

    Member functions are provided to solve linear systems using the LLT
    decomposition.

SourceFiles
    LLTMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef LLTMatrix_H
#define LLTMatrix_H

#include "SquareMatrix.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LLTMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LLTMatrix
:
    public SquareMatrix<Type>
{

public:

    // Constructors

        //- Construct null
        LLTMatrix();

        //- Construct from a square matrix and perform the decomposition
        LLTMatrix(const SquareMatrix<Type>& M);


    // Member Functions

        //- Perform the Cholesky decomposition of the matrix
        void decompose(const SquareMatrix<Type>& M);

        //- Solve the linear system with the given source
        //  and returning the solution in the Field argument x.
        //  This function may be called with the same field for x and source.
        void solve(Field<Type>& x, const Field<Type>& source) const;

        //- Solve the linear system with the given source
        //  returning the solution
        tmp<Field<Type>> solve(const Field<Type>& source) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LLTMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
