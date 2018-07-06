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

Class
    Foam::scalarMatrices

Description
    Scalar matrices

    LUDecompose for scalarSymmetricSquareMatrix implements the Cholesky
    decomposition method from JAMA, a public-domain library developed at NIST,
    available at http://math.nist.gov/tnt/index.html

SourceFiles
    scalarMatrices.C
    scalarMatricesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef scalarMatrices_H
#define scalarMatrices_H

#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "SymmetricSquareMatrix.H"
#include "DiagonalMatrix.H"
#include "scalarField.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef RectangularMatrix<scalar> scalarRectangularMatrix;
typedef SquareMatrix<scalar> scalarSquareMatrix;
typedef SymmetricSquareMatrix<scalar> scalarSymmetricSquareMatrix;
typedef DiagonalMatrix<scalar> scalarDiagonalMatrix;

//- Solve the matrix using Gaussian elimination with pivoting,
//  returning the solution in the source
template<class Type>
void solve(scalarSquareMatrix& matrix, List<Type>& source);

//- Solve the matrix using Gaussian elimination with pivoting
//  and return the solution
template<class Type>
void solve
(
    List<Type>& psi,
    const scalarSquareMatrix& matrix,
    const List<Type>& source
);

//- LU decompose the matrix with pivoting
void LUDecompose
(
    scalarSquareMatrix& matrix,
    labelList& pivotIndices
);

//- LU decompose the matrix with pivoting.
//  sign is -1 for odd number of row interchanges and 1 for even number.
void LUDecompose
(
    scalarSquareMatrix& matrix,
    labelList& pivotIndices,
    label& sign
);

//- LU decompose the matrix into a lower (L) and upper (U) part. U = L.T()
void LUDecompose(scalarSymmetricSquareMatrix& matrix);

//- LU back-substitution with given source, returning the solution
//  in the source
template<class Type>
void LUBacksubstitute
(
    const scalarSquareMatrix& luMmatrix,
    const labelList& pivotIndices,
    List<Type>& source
);

//- LU back-substitution with given source, returning the solution
//  in the source. Specialised for symmetric square matrices that have been
//  decomposed into LU where U = L.T() as pivoting is not required
template<class Type>
void LUBacksubstitute
(
    const scalarSymmetricSquareMatrix& luMmatrix,
    List<Type>& source
);

//- Solve the matrix using LU decomposition with pivoting
//  returning the LU form of the matrix and the solution in the source
template<class Type>
void LUsolve(scalarSquareMatrix& matrix, List<Type>& source);

//- Solve the matrix using LU decomposition returning the LU form of the matrix
//  and the solution in the source, where U = L.T()
template<class Type>
void LUsolve(scalarSymmetricSquareMatrix& matrix, List<Type>& source);

template<class Form, class Type>
void multiply
(
    Matrix<Form, Type>& answer,         // value changed in return
    const Matrix<Form, Type>& A,
    const Matrix<Form, Type>& B
);

void multiply
(
    scalarRectangularMatrix& answer,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B,
    const scalarRectangularMatrix& C
);

void multiply
(
    scalarRectangularMatrix& answer,         // value changed in return
    const scalarRectangularMatrix& A,
    const DiagonalMatrix<scalar>& B,
    const scalarRectangularMatrix& C
);

void multiply
(
    scalarSquareMatrix& answer,         // value changed in return
    const scalarSquareMatrix& A,
    const DiagonalMatrix<scalar>& B,
    const scalarSquareMatrix& C
);

//- Return the inverse of matrix A using SVD
scalarRectangularMatrix SVDinv
(
    const scalarRectangularMatrix& A,
    scalar minCondition = 0
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "scalarMatricesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
