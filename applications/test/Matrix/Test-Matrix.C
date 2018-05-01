/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "scalarMatrices.H"
#include "LUscalarMatrix.H"
#include "LLTMatrix.H"
#include "QRMatrix.H"
#include "vector.H"
#include "tensor.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    SquareMatrix<scalar> hmm(3);

    hmm(0, 0) = -3.0;
    hmm(0, 1) = 10.0;
    hmm(0, 2) = -4.0;
    hmm(1, 0) = 2.0;
    hmm(1, 1) = 3.0;
    hmm(1, 2) = 10.0;
    hmm(2, 0) = 2.0;
    hmm(2, 1) = 6.0;
    hmm(2, 2) = 1.0;

    // Info<< hmm << endl << hmm - 2.0*(-hmm) << endl;
    Info<< max(hmm) << endl;
    Info<< min(hmm) << endl;

    SquareMatrix<scalar> hmm2(3, I);

    hmm = hmm2;

    Info<< hmm << endl;

    // SquareMatrix<scalar> hmm3(Sin);

    // Info<< hmm3 << endl;

    SquareMatrix<scalar> hmm4;

    hmm4 = hmm2;

    Info<< hmm4 << endl;

    SquareMatrix<scalar> hmm5;

    hmm4 = hmm5;
    Info<< hmm5 << endl;

    {
        RectangularMatrix<scalar> rm1(5, 6, 3.1);
        rm1(0, 1) = 4.5;
        RectangularMatrix<scalar> rm1b(rm1.block(2, 2, 0, 0));
        Info<< "rm1b = " << rm1b << endl;
    }

    {
        scalarSymmetricSquareMatrix symmMatrix(3, Zero);

        symmMatrix(0, 0) = 4;
        symmMatrix(1, 0) = 12;
        symmMatrix(1, 1) = 37;
        symmMatrix(2, 0) = -16;
        symmMatrix(2, 1) = -43;
        symmMatrix(2, 2) = 98;

        Info<< "Symmetric Square Matrix = " << symmMatrix << endl;

        Info<< "Inverse = " << inv(symmMatrix) << endl;
        Info<< "Determinant = " << det(symmMatrix) << endl;

        scalarSymmetricSquareMatrix symmMatrix2(symmMatrix);
        LUDecompose(symmMatrix2);

        Info<< "Inverse = " << invDecomposed(symmMatrix2) << endl;
        Info<< "Determinant = " << detDecomposed(symmMatrix2) << endl;

        scalarDiagonalMatrix rhs(3, 0);
        rhs[0] = 1;
        rhs[1] = 2;
        rhs[2] = 3;

        LUsolve(symmMatrix, rhs);

        Info<< "Decomposition = " << symmMatrix << endl;
        Info<< "Solution = " << rhs << endl;
    }


    scalarSquareMatrix squareMatrix(3, Zero);

    squareMatrix(0, 0) = 4;
    squareMatrix(0, 1) = 12;
    squareMatrix(0, 2) = -16;
    squareMatrix(1, 0) = 12;
    squareMatrix(1, 1) = 37;
    squareMatrix(1, 2) = -43;
    squareMatrix(2, 0) = -16;
    squareMatrix(2, 1) = -43;
    squareMatrix(2, 2) = 98;

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

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
