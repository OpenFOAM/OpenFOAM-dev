/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    const label n = matrix.n();

    SymmetricSquareMatrix<Type> inv(n, Zero);

    for (label i=0; i<n; i++)
    {
        inv(i, i) = 1.0/matrix(i, i);

        for (label j=0; j<i; j++)
        {
            Type sum = Zero;

            for (label k=j; k<i; k++)
            {
                sum -= matrix(i, k)*inv(k, j);
            }

            inv(i, j) = sum/matrix(i, i);
        }
    }

    SymmetricSquareMatrix<Type> result(n, Zero);

    for (label k=0; k<n; k++)
    {
        for (label i=0; i <= k; i++)
        {
            for (label j=0; j <= k; j++)
            {
                result(i, j) += inv(k, i)*inv(k, j);
            }
        }
    }

    return result;
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
Type Foam::detDecomposed(const SymmetricSquareMatrix<Type>& matrix)
{
    Type diagProduct = pTraits<Type>::one;

    for (label i=0; i<matrix.m(); i++)
    {
        diagProduct *= matrix(i, i);
    }

    return sqr(diagProduct);
}


template<class Type>
Type Foam::det(const SymmetricSquareMatrix<Type>& matrix)
{
    SymmetricSquareMatrix<Type> matrixTmp(matrix);
    LUDecompose(matrixTmp);

    return detDecomposed(matrixTmp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
