/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "scalarMatrices.H"
#include "Swap.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::solve
(
    scalarSquareMatrix& A,
    List<Type>& sourceSol
)
{
    const label m = A.m();

    // Elimination
    for (label i=0; i<m; i++)
    {
        label iMax = i;
        scalar largestCoeff = mag(A(iMax, i));

        // Swap elements around to find a good pivot
        for (label j=i+1; j<m; j++)
        {
            if (mag(A(j, i)) > largestCoeff)
            {
                iMax = j;
                largestCoeff = mag(A(iMax, i));
            }
        }

        if (i != iMax)
        {
            for (label k=i; k<m; k++)
            {
                Swap(A(i, k), A(iMax, k));
            }
            Swap(sourceSol[i], sourceSol[iMax]);
        }

        // Check that the system of equations isn't singular
        if (mag(A(i, i)) < small)
        {
            FatalErrorInFunction << "Singular Matrix" << exit(FatalError);
        }

        // Reduce to upper triangular form
        for (label j=i+1; j<m; j++)
        {
            const scalar multiplier = A(j, i)/A(i, i);
            sourceSol[j] -= multiplier*sourceSol[i];

            for (label k=m-1; k>=i; k--)
            {
                A(j, k) -= multiplier*A(i, k);
            }
        }
    }

    // Back-substitution
    for (label j=m-1; j>=0; j--)
    {
        Type a = Zero;

        for (label k=j+1; k<m; k++)
        {
            a += A(j, k)*sourceSol[k];
        }

        sourceSol[j] = (sourceSol[j] - a)/A(j, j);
    }
}


template<class Type>
void Foam::solve
(
    List<Type>& psi,
    const scalarSquareMatrix& matrix,
    const List<Type>& source
)
{
    scalarSquareMatrix A = matrix;
    psi = source;
    solve(A, psi);
}


template<class Type>
void Foam::LUBacksubstitute
(
    const scalarSquareMatrix& luMatrix,
    const labelList& pivotIndices,
    List<Type>& sourceSol
)
{
    const label m = luMatrix.m();

    label ii = 0;

    for (label i=0; i<m; i++)
    {
        label ip = pivotIndices[i];
        Type sum = sourceSol[ip];
        sourceSol[ip] = sourceSol[i];
        const scalar* __restrict__ luMatrixi = luMatrix[i];

        if (ii != 0)
        {
            for (label j=ii-1; j<i; j++)
            {
                sum -= luMatrixi[j]*sourceSol[j];
            }
        }
        else if (sum != pTraits<Type>::zero)
        {
            ii = i+1;
        }

        sourceSol[i] = sum;
    }

    for (label i=m-1; i>=0; i--)
    {
        Type sum = sourceSol[i];
        const scalar* __restrict__ luMatrixi = luMatrix[i];

        for (label j=i+1; j<m; j++)
        {
            sum -= luMatrixi[j]*sourceSol[j];
        }

        sourceSol[i] = sum/luMatrixi[i];
    }
}


template<class Type>
void Foam::LUBacksubstitute
(
    const scalarSquareMatrix& luMatrix,
    List<Type>& sourceSol
)
{
    const label m = luMatrix.m();

    label ii = 0;

    for (label i=0; i<m; i++)
    {
        Type sum = sourceSol[i];

        if (ii != 0)
        {
            for (label j=ii-1; j<i; j++)
            {
                sum -= luMatrix(i, j)*sourceSol[j];
            }
        }
        else if (sum != pTraits<Type>::zero)
        {
            ii = i+1;
        }

        sourceSol[i] = sum;
    }

    for (label i=m-1; i>=0; i--)
    {
        Type sum = sourceSol[i];

        for (label j=i+1; j<m; j++)
        {
            sum -= luMatrix(i, j)*sourceSol[j];
        }

        sourceSol[i] = sum/luMatrix(i, i);
    }
}


template<class Type>
void Foam::LUBacksubstitute
(
    const scalarSymmetricSquareMatrix& luMatrix,
    List<Type>& sourceSol
)
{
    const label m = luMatrix.m();

    label ii = 0;

    for (label i=0; i<m; i++)
    {
        Type sum = sourceSol[i];

        if (ii != 0)
        {
            for (label j=ii-1; j<i; j++)
            {
                sum -= luMatrix(i, j)*sourceSol[j];
            }
        }
        else if (sum != pTraits<Type>::zero)
        {
            ii = i+1;
        }

        sourceSol[i] = sum/luMatrix(i, i);
    }

    for (label i=m-1; i>=0; i--)
    {
        Type sum = sourceSol[i];

        for (label j=i+1; j<m; j++)
        {
            sum -= luMatrix(i, j)*sourceSol[j];
        }

        sourceSol[i] = sum/luMatrix(i, i);
    }
}


template<class Type>
void Foam::LUsolve
(
    scalarSquareMatrix& matrix,
    List<Type>& sourceSol
)
{
    labelList pivotIndices(matrix.m());
    scalarList scale(matrix.m());
    LUDecompose(matrix, pivotIndices, scale);
    LUBacksubstitute(matrix, pivotIndices, sourceSol);
}


template<class Type>
void Foam::LUsolve
(
    scalarSymmetricSquareMatrix& matrix,
    List<Type>& sourceSol
)
{
    LUDecompose(matrix);
    LUBacksubstitute(matrix, sourceSol);
}


template<class Form, class Type>
void Foam::multiply
(
    Matrix<Form, Type>& ans,         // value changed in return
    const Matrix<Form, Type>& A,
    const Matrix<Form, Type>& B
)
{
    if (A.n() != B.m())
    {
        FatalErrorInFunction
            << "A and B must have identical inner dimensions but A.n = "
            << A.n() << " and B.m = " << B.m()
            << abort(FatalError);
    }

    ans = Matrix<Form, Type>(A.m(), B.n(), scalar(0));

    for (label i=0; i<A.m(); i++)
    {
        for (label j=0; j<B.n(); j++)
        {
            for (label l=0; l<B.m(); l++)
            {
                ans(i, j) += A(i, l)*B(l, j);
            }
        }
    }
}


// ************************************************************************* //
