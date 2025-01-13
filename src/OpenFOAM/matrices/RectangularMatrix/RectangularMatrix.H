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

Class
    Foam::RectangularMatrix

Description
    A templated 2D rectangular m x n matrix of objects of \<Type\>.

    The matrix dimensions are used for subscript bounds checking etc.

SourceFiles
    RectangularMatrixI.H
    RectangularMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef RectangularMatrix_H
#define RectangularMatrix_H

#include "Matrix.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RectangularMatrix
:
    public Matrix<RectangularMatrix<Type>, Type>
{

public:

    // Constructors

        using Matrix<RectangularMatrix<Type>, Type>::Matrix;

        //- Null constructor.
        inline RectangularMatrix();

        //- Construct as copy of a square matrix
        inline RectangularMatrix(const SquareMatrix<Type>&);

        //- Clone
        inline autoPtr<RectangularMatrix<Type>> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return word("RectangularMatrix<") + pTraits<Type>::typeName + '>';
        }


    // Member Operators

        //- Assignment of all elements to zero
        void operator=(const zero);
};


// Global functions and operators

template<class Type>
class typeOfInnerProduct<Type, RectangularMatrix<Type>, RectangularMatrix<Type>>
{
public:

    typedef RectangularMatrix<Type> type;
};

template<class Type>
class typeOfInnerProduct<Type, RectangularMatrix<Type>, SquareMatrix<Type>>
{
public:

    typedef RectangularMatrix<Type> type;
};

template<class Type>
class typeOfInnerProduct<Type, SquareMatrix<Type>, RectangularMatrix<Type>>
{
public:

    typedef RectangularMatrix<Type> type;
};


template<class Type>
RectangularMatrix<Type> outer(const Field<Type>& f1, const Field<Type>& f2);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RectangularMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
