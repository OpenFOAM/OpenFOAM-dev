/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "SymmetricSquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::SymmetricSquareMatrix<Type> Foam::invDecomposed
(
    const SymmetricSquareMatrix<Type>& matrix
)
{
    SymmetricSquareMatrix<Type> inv(matrix.n(), matrix.n(), 0.0);

    for (label i = 0; i < matrix.n(); ++i)
    {
        inv[i][i] = 1.0/matrix[i][i];

        for (label j = 0; j < i; ++j)
        {
            scalar sum = 0.0;

            for (label k = j; k < i; k++)
            {
                sum -= matrix[i][k]*inv[k][j];
            }

            inv[i][j] = sum/matrix[i][i];
        }
    }

    return inv.T()*inv;
}


template<class Type>
Foam::SymmetricSquareMatrix<Type> Foam::inv
(
    const SymmetricSquareMatrix<Type>& matrix
)
{
    SymmetricSquareMatrix<Type> matrixTmp(matrix);

    LUDecompose(matrixTmp);

    return invDecomposed(matrixTmp);
}


template<class Type>
Foam::scalar Foam::detDecomposed(const SymmetricSquareMatrix<Type>& matrix)
{
    scalar diagProduct = 1.0;

    for (label i = 0; i < matrix.n(); ++i)
    {
        diagProduct *= matrix[i][i];
    }

    return sqr(diagProduct);
}


template<class Type>
Foam::scalar Foam::det(const SymmetricSquareMatrix<Type>& matrix)
{
    SymmetricSquareMatrix<Type> matrixTmp = matrix;

    LUDecompose(matrixTmp);

    return detDecomposed(matrixTmp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
