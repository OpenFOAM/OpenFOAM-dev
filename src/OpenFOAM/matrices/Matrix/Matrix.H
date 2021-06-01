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
    Foam::Matrix

Description
    A templated (m x n) matrix of objects of \<T\>.

SourceFiles
    Matrix.C
    MatrixI.H
    MatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef Matrix_H
#define Matrix_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "Field.H"
#include "zero.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Form, class Type> class Matrix;

template<class Form, class Type> Istream& operator>>
(
    Istream&,
    Matrix<Form, Type>&
);

template<class Form, class Type> Ostream& operator<<
(
    Ostream&,
    const Matrix<Form, Type>&
);

template<class MatrixType>
class ConstMatrixBlock;

template<class MatrixType>
class MatrixBlock;


/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Type>
class Matrix
{
    // Private Data

        //- Number of rows and columns in Matrix.
        label mRows_, nCols_;

        //- Row pointers
        Type* __restrict__ v_;

        //- Allocate the storage for the element vector
        void allocate();


public:

    //- Matrix type
    typedef Matrix<Form, Type> mType;

    //- Component type
    typedef Type cmptType;


    // Static Member Functions

        //- Return a null Matrix
        inline static const mType& null();


    // Constructors

        //- Null constructor.
        inline Matrix();

        //- Construct given number of rows and columns.
        Matrix(const label m, const label n);

        //- Construct with given number of rows and columns
        //  initialising all elements to zero
        Matrix(const label m, const label n, const zero);

        //- Construct with given number of rows and columns
        //  initialising all elements to the given value
        Matrix(const label m, const label n, const Type&);

        //- Copy constructor.
        Matrix(const mType&);

        //- Copy constructor from matrix of a different form
        template<class Form2>
        explicit Matrix(const Matrix<Form2, Type>&);

        //- Construct from a block of another matrix
        template<class MatrixType>
        Matrix(const ConstMatrixBlock<MatrixType>&);

        //- Construct from a block of another matrix
        template<class MatrixType>
        Matrix(const MatrixBlock<MatrixType>&);

        //- Construct from Istream.
        Matrix(Istream&);

        //- Clone
        inline autoPtr<mType> clone() const;


    //- Destructor
    ~Matrix();


    // Member Functions

        // Access

            //- Return the number of rows
            inline label m() const;

            //- Return the number of columns
            inline label n() const;

            //- Return the number of elements in matrix (m*n)
            inline label size() const;

            //- Return element vector of the constant Matrix
            inline const Type* v() const;

            //- Return element vector of the Matrix
            inline Type* v();


        // Block access

            inline ConstMatrixBlock<mType> block
            (
                const label m,
                const label n,
                const label mStart,
                const label nStart
            ) const;

            template<class VectorSpace>
            inline ConstMatrixBlock<mType> block
            (
                const label mStart,
                const label nStart
            ) const;

            inline ConstMatrixBlock<mType> col
            (
                const label m,
                const label rowStart
            ) const;

            inline ConstMatrixBlock<mType> col
            (
                const label m,
                const label mStart,
                const label nStart
            ) const;


            inline MatrixBlock<mType> block
            (
                const label m,
                const label n,
                const label mStart,
                const label nStart
            );

            template<class VectorSpace>
            inline MatrixBlock<mType> block
            (
                const label mStart,
                const label nStart
            );

            inline MatrixBlock<mType> col
            (
                const label m,
                const label rowStart
            );

            inline MatrixBlock<mType> col
            (
                const label m,
                const label mStart,
                const label nStart
            );


        // Check

            //- Check index i is within valid range (0 ... m-1).
            inline void checki(const label i) const;

            //- Check index j is within valid range (0 ... n-1).
            inline void checkj(const label j) const;


        // Edit

            //- Clear the Matrix, i.e. set sizes to zero.
            void clear();

            //- Transfer the contents of the argument Matrix into this Matrix
            //  and annul the argument Matrix.
            void transfer(mType&);

            //- Resize the matrix preserving the elements
            void setSize(const label m, const label n);

            //- Resize the matrix without reallocating storage (unsafe)
            inline void shallowResize(const label m, const label n);


        //- Return the transpose of the matrix
        Form T() const;


    // Member Operators

        //- Return subscript-checked row of Matrix.
        inline Type* operator[](const label);

        //- Return subscript-checked row of constant Matrix.
        inline const Type* operator[](const label) const;

        //- (i, j) const element access operator
        inline const Type& operator()(const label i, const label j) const;

        //- (i, j) element access operator
        inline Type& operator()(const label i, const label j);

        //- Assignment operator. Takes linear time.
        void operator=(const mType&);

        //- Assignment to a block of another matrix
        template<class MatrixType>
        void operator=(const ConstMatrixBlock<MatrixType>&);

        //- Assignment to a block of another matrix
        template<class MatrixType>
        void operator=(const MatrixBlock<MatrixType>&);

        //- Assignment of all elements to zero
        void operator=(const zero);

        //- Assignment of all elements to the given value
        void operator=(const Type&);


    // IOstream Operators

        //- Read Matrix from Istream, discarding contents of existing Matrix.
        friend Istream& operator>> <Form, Type>
        (
            Istream&,
            mType&
        );

        //- Write Matrix to Ostream.
        friend Ostream& operator<< <Form, Type>
        (
            Ostream&,
            const mType&
        );
};


// Global functions and operators

template<class Form, class Type>
void writeEntry(Ostream&, const Matrix<Form, Type>&);

template<class Form, class Type>
const Type& max(const Matrix<Form, Type>&);

template<class Form, class Type>
const Type& min(const Matrix<Form, Type>&);

template<class Form, class Type>
Form operator-(const Matrix<Form, Type>&);

template<class Form, class Type>
Form operator+
(
    const Matrix<Form, Type>&,
    const Matrix<Form, Type>&
);

template<class Form, class Type>
Form operator-
(
    const Matrix<Form, Type>&,
    const Matrix<Form, Type>&
);

template<class Form, class Type>
Form operator*
(
    const scalar,
    const Matrix<Form, Type>&
);

template<class Form, class Type>
Form operator*
(
    const Matrix<Form, Type>&,
    const scalar
);

template<class Form, class Type>
Form operator/
(
    const Matrix<Form, Type>&,
    const scalar
);

template<class Form1, class Form2, class Type>
typename typeOfInnerProduct<Type, Form1, Form2>::type
operator*
(
    const Matrix<Form1, Type>& a,
    const Matrix<Form2, Type>& b
);

template<class Form, class Type>
tmp<Field<Type>> operator*
(
    const Matrix<Form, Type>&,
    const Field<Type>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Matrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
