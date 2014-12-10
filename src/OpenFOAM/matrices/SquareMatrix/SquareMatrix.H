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

Class
    Foam::SquareMatrix

Description
    A templated 2D square matrix of objects of \<T\>, where the n x n matrix
    dimension is known and used for subscript bounds checking, etc.

SourceFiles
    SquareMatrixI.H
    SquareMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef SquareMatrix_H
#define SquareMatrix_H

#include "Matrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SquareMatrix
:
    public Matrix<SquareMatrix<Type>, Type>
{

public:

    // Constructors

        //- Null constructor.
        inline SquareMatrix();

        //- Construct given number of rows/columns.
        inline SquareMatrix(const label n);

        //- Construct given number of rows and columns,
        //  It checks that m == n.
        inline SquareMatrix(const label m, const label n);

        //- Construct with given number of rows and rows
        //  and value for all elements.
        //  It checks that m == n.
        inline SquareMatrix(const label m, const label n, const Type&);

        //- Construct from Istream.
        inline SquareMatrix(Istream&);

        //- Clone
        inline autoPtr<SquareMatrix<Type> > clone() const;
};


// Global functions

//- Return the LU decomposed SquareMatrix det
template<class Type>
scalar detDecomposed(const SquareMatrix<Type>&, const label sign);

//- Return the SquareMatrix det
template<class Type>
scalar det(const SquareMatrix<Type>&);

//- Return the SquareMatrix det and the LU decomposition in the original matrix
template<class Type>
scalar det(SquareMatrix<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "SquareMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SquareMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
