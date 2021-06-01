/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "Identity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class RectangularMatrix;


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

        //- Construct from a block of another matrix
        template<class MatrixType>
        inline SquareMatrix(const ConstMatrixBlock<MatrixType>&);

        //- Construct from a block of another matrix
        template<class MatrixType>
        inline SquareMatrix(const MatrixBlock<MatrixType>&);

        //- Construct given number of rows/columns
        //  initialising all elements to zero
        inline SquareMatrix(const label n, const zero);

        //- Construct given number of rows and columns (checked to be equal)
        //  initialising all elements to zero
        inline SquareMatrix(const label m, const label n, const zero);

        //- Construct given number of rows/columns
        //  Initialising to the identity matrix
        inline SquareMatrix(const label n, const Identity<Type>);

        //- Construct with given number of rows and rows
        //  initialising all elements to the given value
        inline SquareMatrix(const label n, const Type&);

        //- Construct as copy of a RectangularMatrix
        //  which is checked to be square
        inline explicit SquareMatrix(const RectangularMatrix<Type>&);

        //- Construct from Istream.
        inline SquareMatrix(Istream&);

        //- Clone
        inline autoPtr<SquareMatrix<Type>> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return word("SquareMatrix<") + pTraits<Type>::typeName + '>';
        }

        // Edit

            //- Resize the matrix preserving the elements
            inline void setSize(const label m);

            //- Resize the matrix without reallocating storage (unsafe)
            inline void shallowResize(const label m);


    // Member Operators

        //- Assignment of all elements to zero
        void operator=(const zero);

        //- Assignment elements to the
        void operator=(const Identity<Type>);
};


// Global functions and operators

//- Return the LU decomposed SquareMatrix det
template<class Type>
scalar detDecomposed(const SquareMatrix<Type>&, const label sign);

//- Return the SquareMatrix det
template<class Type>
scalar det(const SquareMatrix<Type>&);

//- Return the SquareMatrix det and the LU decomposition in the original matrix
template<class Type>
scalar det(SquareMatrix<Type>&);

template<class Type>
class typeOfInnerProduct<Type, SquareMatrix<Type>, SquareMatrix<Type>>
{
public:

    typedef SquareMatrix<Type> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SquareMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SquareMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
