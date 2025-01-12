/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
#include "LUscalarMatrix.H"
#include "LLTMatrix.H"
#include "QRMatrix.H"
#include "eigendecomposition.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    {
        SquareMatrix<scalar> hmm
        {
            {-3.0, 10.0, -4.0},
            {2.0, 3.0, 10.0},
            {2.0, 6.0, 1.0}
        };

        Info<< max(hmm) << endl;
        Info<< min(hmm) << endl;

        SquareMatrix<scalar> hmm2(3, I);

        hmm = hmm2;

        Info<< hmm << endl;

        SquareMatrix<scalar> hmm4;

        hmm4 = hmm2;

        Info<< hmm4 << endl;

        SquareMatrix<scalar> hmm5;

        hmm4 = hmm5;
        Info<< hmm5 << endl;
    }

    {
        RectangularMatrix<scalar> rm1(5, 6, 3.1);
        rm1(0, 1) = 4.5;
        RectangularMatrix<scalar> rm1b(rm1.block(2, 2, 0, 0));
        Info<< "rm1b = " << rm1b << endl;
    }

    {
        scalarSymmetricSquareMatrix symmMatrix
        {
            4, 0, 0, 12, 37, 0, -16, -43, 98
        };

        Info<< "Symmetric Square Matrix = " << symmMatrix << endl;

        Info<< "Inverse = " << inv(symmMatrix) << endl;
        Info<< "Determinant = " << det(symmMatrix) << endl;

        scalarSymmetricSquareMatrix symmMatrix2(symmMatrix);
        LUDecompose(symmMatrix2);

        Info<< "Inverse = " << invDecomposed(symmMatrix2) << endl;
        Info<< "Determinant = " << detDecomposed(symmMatrix2) << endl;

        scalarDiagonalMatrix rhs{1, 2, 3};

        LUsolve(symmMatrix, rhs);

        Info<< "Decomposition = " << symmMatrix << endl;
        Info<< "Solution = " << rhs << endl;
    }


    scalarSquareMatrix squareMatrix{4, 12, -16, 12, 37, -43, -16, -43, 98};

    Info<< nl << "Square Matrix = " << squareMatrix << endl;

    const scalarField source(3, 1);

    {
        {
            scalarSquareMatrix sm(squareMatrix);
            Info<< "det = " << det(sm) << endl;
        }

        scalarSquareMatrix sm(squareMatrix);
        labelList rhs(3, 0);
        label sign;
        LUDecompose(sm, rhs, sign);

        Info<< "Decomposition = " << sm << endl;
        Info<< "Pivots = " << rhs << endl;
        Info<< "Sign = " << sign << endl;
        Info<< "det = " << detDecomposed(sm, sign) << endl;
    }

    {
        LUscalarMatrix LU(squareMatrix);
        scalarField x(LU.solve(source));
        Info<< "LU solve residual " << (squareMatrix*x - source) << endl;

        scalarSquareMatrix inv(3);
        LU.inv(inv);
        Info<< "LU inv " << inv << endl;
        Info<< "LU inv*squareMatrix " << (inv*squareMatrix) << endl;
    }

    {
        LLTMatrix<scalar> LLT(squareMatrix);
        scalarField x(LLT.solve(source));
        Info<< "LLT solve residual " << (squareMatrix*x - source) << endl;
    }

    {
        QRMatrix<scalarSquareMatrix> QR(squareMatrix);
        scalarField x(QR.solve(source));

        Info<< "QR solve residual "
            << (squareMatrix*x - source) << endl;

        Info<< "QR inverse solve residual "
            << (x - QR.inv()*source) << endl;

        Info<< "QR inv *squareMatrix " << (QR.inv()*squareMatrix) << endl;
    }

    {
        scalarSquareMatrix A{1, 0, -4, 0, 5, 4, -4, 4, 3};

        eigendecomposition eigen(A);

        Info<< "matrix " << A << endl;
        Info<< "eigenvalues " << eigen.d() << endl;
        Info<< "eigenvectors " << eigen.V().T() << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
