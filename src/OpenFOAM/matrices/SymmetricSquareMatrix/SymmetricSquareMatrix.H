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
    Foam::SymmetricSquareMatrix

Description
    A templated 2D square symmetric matrix of objects of \<T\>, where the
    n x n matrix dimension is known and used for subscript bounds checking, etc.

SourceFiles
    SymmetricSquareMatrixI.H
    SymmetricSquareMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef SymmetricSquareMatrix_H
#define SymmetricSquareMatrix_H

#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SymmetricSquareMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SymmetricSquareMatrix
:
    public Matrix<SymmetricSquareMatrix<Type>, Type>
{

public:

    // Constructors

        //- Null constructor.
        inline SymmetricSquareMatrix();

        //- Construct given number of rows/columns.
        inline SymmetricSquareMatrix(const label n);

        //- Construct with given number of rows/columns
        inline SymmetricSquareMatrix(const label m, const label n);

        //- Construct with given number of rows/columns
        //  and value for all elements.
        inline SymmetricSquareMatrix(const label m, const label n, const Type&);

        //- Construct from Istream.
        inline SymmetricSquareMatrix(Istream&);

        //- Clone
        inline autoPtr<SymmetricSquareMatrix<Type> > clone() const;


        //- Return subscript-checked row of Matrix.
        inline Type& operator()(const label r, const label c);

        //- Return subscript-checked row of constant Matrix.
        inline const Type& operator()(const label r, const label c) const;
};


// Global functions

//- Return the LU decomposed SymmetricSquareMatrix inverse
template<class Type>
SymmetricSquareMatrix<Type> invDecomposed(const SymmetricSquareMatrix<Type>&);

//- Return the SymmetricSquareMatrix inverse
template<class Type>
SymmetricSquareMatrix<Type> inv(const SymmetricSquareMatrix<Type>&);

//- Return the LU decomposed SymmetricSquareMatrix det
template<class Type>
scalar detDecomposed(const SymmetricSquareMatrix<Type>&);

//- Return the SymmetricSquareMatrix det
template<class Type>
scalar det(const SymmetricSquareMatrix<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "SymmetricSquareMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SymmetricSquareMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
