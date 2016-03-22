/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    A templated 2D matrix of objects of \<T\>, where the n x m matrix
    dimensions are known and used for subscript bounds checking, etc.

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
#include "MatrixSpace.H"
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


/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Type>
class Matrix
{
    // Private data

        //- Number of rows and columns in Matrix.
        label nRows_, nCols_;

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


    // Sub-Block Classes

        class ConstBlock
        {
            // Private data

                //- Const reference to the parent matrix
                const mType& matrix_;

                // Block size
                const label nRows_;
                const label nCols_;

                // Block location in parent matrix
                const label rowStart_;
                const label colStart_;

        public:

            typedef typename mType::cmptType cmptType;

            // Constructors

                //- Construct block for matrix, size and location
                inline ConstBlock
                (
                    const mType& matrix,
                    const label m,
                    const label n,
                    const label mStart,
                    const label nStart
                );


            // Member Functions

                //- Return the number of rows in the block
                inline label m() const;

                //- Return the number of columns in the block
                inline label n() const;

                //- (i, j) const element access operator
                inline const Type& operator()
                (
                    const label i,
                    const label j
                ) const;

                //- Convert a column of a matrix to a Field
                operator Field<Type>() const;
        };


        class Block
        {
            // Private data

                //- Reference to the parent matrix
                mType& matrix_;

                // Block size
                const label nRows_;
                const label nCols_;

                // Block location in parent matrix
                const label rowStart_;
                const label colStart_;

        public:

            typedef typename mType::cmptType cmptType;

            // Constructors

                //- Construct block for matrix, size and location
                inline Block
                (
                    mType& matrix,
                    const label m,
                    const label n,
                    const label mStart,
                    const label nStart
                );


            // Member Functions

                //- Return the number of rows in the block
                inline label m() const;

                //- Return the number of columns in the block
                inline label n() const;

                //- (i, j) const element access operator
                inline const Type& operator()
                (
                    const label i,
                    const label j
                ) const;

                //- (i, j) element access operator
                inline Type& operator()(const label i, const label j);

                //- Convert a column of a matrix to a Field
                operator Field<Type>() const;


            // Member operators

                //- Assignment to a compatible matrix
                void operator=(const mType&);

                //- Assignment to a compatible block
                void operator=(const Block&);

                //- Assignment to a compatible const block
                void operator=(const ConstBlock&);

                //- Assignment to a compatible MatrixSpace
                template<class MSForm, direction Nrows, direction Ncols>
                void operator=(const MatrixSpace<MSForm, Type, Nrows, Ncols>&);

                //- Assignment to a compatible MatrixSpace block
                template
                <
                    template<class, direction, direction> class Block,
                    class SubTensor,
                    direction BRowStart,
                    direction BColStart
                >
                void operator=(const Block<SubTensor, BRowStart, BColStart>&);

                //- Assignment to a compatible VectorSpace (column-vector)
                template<class VSForm, direction Ncmpts>
                void operator=(const VectorSpace<VSForm, Type, Ncmpts>&);

                //- Assignment to a compatible VectorSpace (column-vector) block
                template
                <
                    template<class, direction> class Block,
                    class SubVector,
                    direction BStart
                >
                void operator=(const Block<SubVector, BStart>&);

                //- Assignment to a Field (column-vector)
                void operator=(const Field<Type>&);
        };


    // Constructors

        //- Null constructor.
        inline Matrix();

        //- Construct given number of rows and columns.
        Matrix(const label m, const label n);

        //- Construct with given number of rows and columns
        //  initializing all elements to zero
        Matrix(const label m, const label n, const zero);

        //- Construct with given number of rows and columns
        //  initializing all elements to the given value
        Matrix(const label m, const label n, const Type&);

        //- Copy constructor.
        Matrix(const mType&);

        //- Construct from a block of another matrix
        Matrix(const mType::Block&);

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

            inline ConstBlock block
            (
                const label m,
                const label n,
                const label mStart,
                const label nStart
            ) const;

            template<class VectorSpace>
            inline ConstBlock block
            (
                const label mStart,
                const label nStart
            ) const;

            inline ConstBlock col
            (
                const label m,
                const label rowStart
            ) const;

            inline ConstBlock col
            (
                const label m,
                const label mStart,
                const label nStart
            ) const;


            inline Block block
            (
                const label m,
                const label n,
                const label mStart,
                const label nStart
            );

            template<class VectorSpace>
            inline Block block(const label mStart, const label nStart);

            inline Block col(const label m, const label rowStart);

            inline Block col
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
            void setSize(const label m, const label nCols);


        //- Return the transpose of the matrix
        Form T() const;


    // Member operators

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
        void operator=(const mType::Block&);

        //- Assignment of all entries to zero
        void operator=(const zero);

        //- Assignment of all entries to the given value
        void operator=(const Type&);


    // IOstream operators

        //- Read Matrix from Istream, discarding contents of existing Matrix.
        friend Istream& operator>> <Form, Type>
        (
            Istream&,
            mType&
        );

        // Write Matrix to Ostream.
        friend Ostream& operator<< <Form, Type>
        (
            Ostream&,
            const mType&
        );
};


// Global functions and operators

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

template<class Form, class Type>
Form operator*
(
    const Matrix<Form, Type>&,
    const Matrix<Form, Type>&
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
