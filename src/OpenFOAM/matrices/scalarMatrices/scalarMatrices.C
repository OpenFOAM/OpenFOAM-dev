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

\*---------------------------------------------------------------------------*/

#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineCompoundTypeName(RectangularMatrix<scalar>, scalarRectangularMatrix);
addCompoundToRunTimeSelectionTable
(
    RectangularMatrix<scalar>,
    scalarRectangularMatrix
);

defineCompoundTypeName(SquareMatrix<scalar>, scalarSquareMatrix);
addCompoundToRunTimeSelectionTable(SquareMatrix<scalar>, scalarSquareMatrix);

defineCompoundTypeName
(
    SymmetricSquareMatrix<scalar>,
    scalarSymmetricSquareMatrix
);
addCompoundToRunTimeSelectionTable
(
    SymmetricSquareMatrix<scalar>,
    scalarSymmetricSquareMatrix
);

defineCompoundTypeName(DiagonalMatrix<scalar>, scalarDiagonalMatrix);
addCompoundToRunTimeSelectionTable
(
    DiagonalMatrix<scalar>,
    scalarDiagonalMatrix
);

}

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
    label m = matrix.m();
    scalar vv[m];
    sign = 1;

    for (label i=0; i<m; i++)
    {
        scalar largestCoeff = 0.0;
        scalar temp;
        const scalar* __restrict__ matrixi = matrix[i];

        for (label j=0; j<m; j++)
        {
            if ((temp = mag(matrixi[j])) > largestCoeff)
            {
                largestCoeff = temp;
            }
        }

        if (largestCoeff == 0.0)
        {
            FatalErrorInFunction
                << "Singular matrix" << exit(FatalError);
        }

        vv[i] = 1.0/largestCoeff;
    }

    for (label j=0; j<m; j++)
    {
        scalar* __restrict__ matrixj = matrix[j];

        for (label i=0; i<j; i++)
        {
            scalar* __restrict__ matrixi = matrix[i];

            scalar sum = matrixi[j];
            for (label k=0; k<i; k++)
            {
                sum -= matrixi[k]*matrix(k, j);
            }
            matrixi[j] = sum;
        }

        label iMax = 0;

        scalar largestCoeff = 0.0;
        for (label i=j; i<m; i++)
        {
            scalar* __restrict__ matrixi = matrix[i];
            scalar sum = matrixi[j];

            for (label k=0; k<j; k++)
            {
                sum -= matrixi[k]*matrix(k, j);
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

            for (label k=0; k<m; k++)
            {
                Swap(matrixj[k], matrixiMax[k]);
            }

            sign *= -1;
            vv[iMax] = vv[j];
        }

        if (matrixj[j] == 0.0)
        {
            matrixj[j] = small;
        }

        if (j != m-1)
        {
            scalar rDiag = 1.0/matrixj[j];

            for (label i=j+1; i<m; i++)
            {
                matrix(i, j) *= rDiag;
            }
        }
    }
}


void Foam::LUDecompose(scalarSymmetricSquareMatrix& matrix)
{
    // Store result in upper triangular part of matrix
    label size = matrix.m();

    // Set upper triangular parts to zero.
    for (label j=0; j<size; j++)
    {
        for (label k=j + 1; k<size; k++)
        {
            matrix(j, k) = 0.0;
        }
    }

    for (label j=0; j<size; j++)
    {
        scalar d = 0.0;

        for (label k=0; k<j; k++)
        {
            scalar s = 0.0;

            for (label i=0; i<k; i++)
            {
                s += matrix(i, k)*matrix(i, j);
            }

            s = (matrix(j, k) - s)/matrix(k, k);

            matrix(k, j) = s;
            matrix(j, k) = s;

            d += sqr(s);
        }

        d = matrix(j, j) - d;

        if (d < 0.0)
        {
            FatalErrorInFunction
                << "Matrix is not symmetric positive-definite. Unable to "
                << "decompose."
                << abort(FatalError);
        }

        matrix(j, j) = sqrt(d);
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
    if (A.n() != B.m())
    {
        FatalErrorInFunction
            << "A and B must have identical inner dimensions but A.n = "
            << A.n() << " and B.m = " << B.m()
            << abort(FatalError);
    }

    if (B.n() != C.m())
    {
        FatalErrorInFunction
            << "B and C must have identical inner dimensions but B.n = "
            << B.n() << " and C.m = " << C.m()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.m(), C.n(), scalar(0));

    for (label i=0; i<A.m(); i++)
    {
        for (label g = 0; g < C.n(); g++)
        {
            for (label l=0; l<C.m(); l++)
            {
                scalar ab = 0;
                for (label j=0; j<A.n(); j++)
                {
                    ab += A(i, j)*B(j, l);
                }
                ans(i, g) += C(l, g) * ab;
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
    if (A.n() != B.size())
    {
        FatalErrorInFunction
            << "A and B must have identical inner dimensions but A.n = "
            << A.n() << " and B.m = " << B.size()
            << abort(FatalError);
    }

    if (B.size() != C.m())
    {
        FatalErrorInFunction
            << "B and C must have identical inner dimensions but B.n = "
            << B.size() << " and C.m = " << C.m()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.m(), C.n(), scalar(0));

    for (label i=0; i<A.m(); i++)
    {
        for (label g=0; g<C.n(); g++)
        {
            for (label l=0; l<C.m(); l++)
            {
                ans(i, g) += C(l, g) * A(i, l)*B[l];
            }
        }
    }
}


void Foam::multiply
(
    scalarSquareMatrix& ans,         // value changed in return
    const scalarSquareMatrix& A,
    const DiagonalMatrix<scalar>& B,
    const scalarSquareMatrix& C
)
{
    if (A.m() != B.size())
    {
        FatalErrorInFunction
            << "A and B must have identical dimensions but A.m = "
            << A.m() << " and B.m = " << B.size()
            << abort(FatalError);
    }

    if (B.size() != C.m())
    {
        FatalErrorInFunction
            << "B and C must have identical dimensions but B.m = "
            << B.size() << " and C.m = " << C.m()
            << abort(FatalError);
    }

    const label size = A.m();

    ans = scalarSquareMatrix(size, Zero);

    for (label i=0; i<size; i++)
    {
        for (label g=0; g<size; g++)
        {
            for (label l=0; l<size; l++)
            {
                ans(i, g) += C(l, g)*A(i, l)*B[l];
            }
        }
    }
}


Foam::scalarRectangularMatrix Foam::SVDinv
(
    const scalarRectangularMatrix& A,
    scalar minCondition
)
{
    SVD svd(A, minCondition);
    return svd.VSinvUt();
}


// ************************************************************************* //
