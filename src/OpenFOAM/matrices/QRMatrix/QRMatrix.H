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
    Foam::QRMatrix

Description
    Class templated on matrix type to perform the QR decomposition using
    Householder reflections on a square or rectangular matrix.

SourceFiles
    QRMatrixI.H
    QRMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef QRMatrix_H
#define QRMatrix_H

#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class QRMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class QRMatrix
{

public:

    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<cmptType> QMatrixType;
    typedef MatrixType RMatrixType;

private:

    // Private data

        //- The Q-matrix
        QMatrixType Q_;

        //- The R-matrix
        RMatrixType R_;


    // Private member functions

        //- Solve the linear system with the Field argument x initialized to
        //  the appropriate transformed source (e.g. Q.T()*source)
        //  and return the solution in x
        void solvex(Field<cmptType>& x) const;


public:

    // Constructors

        //- Construct null
        inline QRMatrix();

        //- Construct decomposing given matrix
        QRMatrix(const MatrixType& M);


    // Member Functions

        //- Return Q-matrix
        inline const QMatrixType& Q() const;

        //- Return R-matrix
        inline const RMatrixType& R() const;

        //- Decompose given matrix
        void decompose(const MatrixType& M);

        //- Solve the linear system with the given source
        //  and returning the solution in the Field argument x
        void solve(Field<cmptType>& x, const Field<cmptType>& source) const;

        //- Solve the linear system with the given source
        //  returning the solution
        tmp<Field<cmptType>> solve(const Field<cmptType>& source) const;

        //- Return the inverse of a square matrix
        QMatrixType inv() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "QRMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "QRMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
