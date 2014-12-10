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

#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LUDecompose
(
    scalarSquareMatrix& matrix,
    labelList& pivotIndices
)
{
    label sign;
    LUDecompose(matrix, pivotIndices, sign);
}


void Foam::LUDecompose
(
    scalarSquareMatrix& matrix,
    labelList& pivotIndices,
    label& sign
)
{
    label n = matrix.n();
    scalar vv[n];
    sign = 1;

    for (register label i=0; i<n; i++)
    {
        scalar largestCoeff = 0.0;
        scalar temp;
        const scalar* __restrict__ matrixi = matrix[i];

        for (register label j=0; j<n; j++)
        {
            if ((temp = mag(matrixi[j])) > largestCoeff)
            {
                largestCoeff = temp;
            }
        }

        if (largestCoeff == 0.0)
        {
            FatalErrorIn
            (
                "LUdecompose"
                "(scalarSquareMatrix& matrix, labelList& rowIndices)"
            )   << "Singular matrix" << exit(FatalError);
        }

        vv[i] = 1.0/largestCoeff;
    }

    for (register label j=0; j<n; j++)
    {
        scalar* __restrict__ matrixj = matrix[j];

        for (register label i=0; i<j; i++)
        {
            scalar* __restrict__ matrixi = matrix[i];

            scalar sum = matrixi[j];
            for (register label k=0; k<i; k++)
            {
                sum -= matrixi[k]*matrix[k][j];
            }
            matrixi[j] = sum;
        }

        label iMax = 0;

        scalar largestCoeff = 0.0;
        for (register label i=j; i<n; i++)
        {
            scalar* __restrict__ matrixi = matrix[i];
            scalar sum = matrixi[j];

            for (register label k=0; k<j; k++)
            {
                sum -= matrixi[k]*matrix[k][j];
            }

            matrixi[j] = sum;

            scalar temp;
            if ((temp = vv[i]*mag(sum)) >= largestCoeff)
            {
                largestCoeff = temp;
                iMax = i;
            }
        }

        pivotIndices[j] = iMax;

        if (j != iMax)
        {
            scalar* __restrict__ matrixiMax = matrix[iMax];

            for (register label k=0; k<n; k++)
            {
                Swap(matrixj[k], matrixiMax[k]);
            }

            sign *= -1;
            vv[iMax] = vv[j];
        }

        if (matrixj[j] == 0.0)
        {
            matrixj[j] = SMALL;
        }

        if (j != n-1)
        {
            scalar rDiag = 1.0/matrixj[j];

            for (register label i=j+1; i<n; i++)
            {
                matrix[i][j] *= rDiag;
            }
        }
    }
}


void Foam::LUDecompose(scalarSymmetricSquareMatrix& matrix)
{
    // Store result in upper triangular part of matrix
    label size = matrix.n();

    // Set upper triangular parts to zero.
    for (label j = 0; j < size; j++)
    {
        for (label k = j + 1; k < size; k++)
        {
            matrix[j][k] = 0.0;
        }
    }

    for (label j = 0; j < size; j++)
    {
        scalar d = 0.0;

        for (label k = 0; k < j; k++)
        {
            scalar s = 0.0;

            for (label i = 0; i < k; i++)
            {
                s += matrix[i][k]*matrix[i][j];
            }

            s = (matrix[j][k] - s)/matrix[k][k];

            matrix[k][j] = s;
            matrix[j][k] = s;

            d += sqr(s);
        }

        d = matrix[j][j] - d;

        if (d < 0.0)
        {
            FatalErrorIn("Foam::LUDecompose(scalarSymmetricSquareMatrix&)")
                << "Matrix is not symmetric positive-definite. Unable to "
                << "decompose."
                << abort(FatalError);
        }

        matrix[j][j] = sqrt(d);
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::multiply
(
    scalarRectangularMatrix& ans,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B,
    const scalarRectangularMatrix& C
)
{
    if (A.m() != B.n())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.m() << " and B.n = " << B.n()
            << abort(FatalError);
    }

    if (B.m() != C.n())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "B and C must have identical inner dimensions but B.m = "
            << B.m() << " and C.n = " << C.n()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.n(), C.m(), scalar(0));

    for (register label i = 0; i < A.n(); i++)
    {
        for (register label g = 0; g < C.m(); g++)
        {
            for (register label l = 0; l < C.n(); l++)
            {
                scalar ab = 0;
                for (register label j = 0; j < A.m(); j++)
                {
                    ab += A[i][j]*B[j][l];
                }
                ans[i][g] += C[l][g] * ab;
            }
        }
    }
}


void Foam::multiply
(
    scalarRectangularMatrix& ans,         // value changed in return
    const scalarRectangularMatrix& A,
    const DiagonalMatrix<scalar>& B,
    const scalarRectangularMatrix& C
)
{
    if (A.m() != B.size())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const DiagonalMatrix<scalar>& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.m() << " and B.n = " << B.size()
            << abort(FatalError);
    }

    if (B.size() != C.n())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const DiagonalMatrix<scalar>& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "B and C must have identical inner dimensions but B.m = "
            << B.size() << " and C.n = " << C.n()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.n(), C.m(), scalar(0));

    for (register label i = 0; i < A.n(); i++)
    {
        for (register label g = 0; g < C.m(); g++)
        {
            for (register label l = 0; l < C.n(); l++)
            {
                ans[i][g] += C[l][g] * A[i][l]*B[l];
            }
        }
    }
}


Foam::RectangularMatrix<Foam::scalar> Foam::SVDinv
(
    const scalarRectangularMatrix& A,
    scalar minCondition
)
{
    SVD svd(A, minCondition);
    return svd.VSinvUt();
}


// ************************************************************************* //
