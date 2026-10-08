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
    scalarSquareMatrix& A,
    labelList& pivotIndices,
    scalarList& scale
)
{
    label sign;
    LUDecompose(A, pivotIndices, scale, sign);
}


void Foam::LUDecompose
(
    scalarSquareMatrix& A,
    labelList& pivotIndices,
    scalarList& scale,
    label& sign
)
{
    const label m = A.m();
    sign = 1;

    // Calculate row scales
    for (label i=0; i<m; i++)
    {
        scalar maxCoeff = 0.0;
        for (label j=0; j<m; j++)
        {
            maxCoeff = max(maxCoeff, mag(A(i, j)));
        }

        if (maxCoeff == 0)
        {
            FatalErrorInFunction << "Singular A" << exit(FatalError);
        }

        scale[i] = maxCoeff;
    }

    // ikj loop with row-swapping
    for (label i=0; i<m; i++)
    {
        // Check in pivoting is required
        label pivotRow = i;
        scalar max_ratio = mag(A(i, i))/scale[i];

        for (label r=i+1; r<m; r++)
        {
            const scalar ratio = mag(A(r, i))/scale[r];
            if (ratio > max_ratio)
            {
                max_ratio = ratio;
                pivotRow = r;
            }
        }

        // Set the pivot row
        pivotIndices[i] = pivotRow;

        // Swap rows if required
        if (pivotRow != i)
        {
            // Swap rows of A
            for (label k=0; k<m; k++)
            {
                Swap(A(i, k), A(pivotRow, k));
            }

            //  Swap the scale elements
            Swap(scale[i], scale[pivotRow]);

            sign *= -1;
        }

        // Check for singularity
        if (mag(A(i, i)) < small)
        {
            FatalErrorInFunction << "Singular A" << exit(FatalError);
        }

        // kj loop
        for (label k=i+1; k<m; k++)
        {
            A(k, i) /= A(i, i);

            for (label j=i+1; j<m; j++)
            {
                A(k, j) -= A(k, i)*A(i, j);
            }
        }
    }
}


void Foam::LUDecompose(scalarSquareMatrix& A, const scalar rowTol)
{
    const label m = A.m();

    for (label i=0; i<m; i++)
    {
        for (label k=0; k<i; k++)
        {
            // Compute the row multiplier for the lower triangular matrix L
            A(i, k) /= A(k, k);

            // If the row multiplier is 0 skip the inner j loop
            if (mag(A(i, k)) <= rowTol)
            {
                A(i, k) = 0;
                continue;
            }

            // Update the remaining elements of the row
            for (label j=k+1; j<m; j++)
            {
                A(i, j) -= A(i, k)*A(k, j);
            }
        }
    }
}


void Foam::LUDecompose(scalarSymmetricSquareMatrix& A)
{
    // Store result in upper triangular part of matrix
    const label m = A.m();

    // Set upper triangular parts to zero.
    for (label j=0; j<m; j++)
    {
        for (label k=j + 1; k<m; k++)
        {
            A(j, k) = 0;
        }
    }

    for (label j=0; j<m; j++)
    {
        scalar d = 0;

        for (label k=0; k<j; k++)
        {
            scalar s = 0;

            for (label i=0; i<k; i++)
            {
                s += A(i, k)*A(i, j);
            }

            s = (A(j, k) - s)/A(k, k);

            A(k, j) = s;
            A(j, k) = s;

            d += sqr(s);
        }

        d = A(j, j) - d;

        if (d < 0)
        {
            FatalErrorInFunction
                << "Matrix is not symmetric positive-definite. Unable to "
                << "decompose."
                << abort(FatalError);
        }

        A(j, j) = sqrt(d);
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
