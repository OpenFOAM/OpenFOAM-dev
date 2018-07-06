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

\*---------------------------------------------------------------------------*/

#include "Matrix.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::allocate()
{
    if (mRows_ && nCols_)
    {
        v_ = new Type[size()];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const zero)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = Zero;
        }
    }
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const Type& s)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = s;
        }
    }
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form, Type>& M)
:
    mRows_(M.mRows_),
    nCols_(M.nCols_),
    v_(nullptr)
{
    if (M.v_)
    {
        allocate();

        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = M.v_[i];
        }
    }
}


template<class Form, class Type>
template<class Form2>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form2, Type>& M)
:
    mRows_(M.m()),
    nCols_(M.n()),
    v_(nullptr)
{
    if (M.v())
    {
        allocate();

        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = M.v()[i];
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
inline Foam::Matrix<Form, Type>::Matrix
(
    const ConstMatrixBlock<MatrixType>& Mb
)
:
    mRows_(Mb.m()),
    nCols_(Mb.n())
{
    allocate();

    for (label i=0; i<mRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i,j) = Mb(i,j);
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
inline Foam::Matrix<Form, Type>::Matrix
(
    const MatrixBlock<MatrixType>& Mb
)
:
    mRows_(Mb.m()),
    nCols_(Mb.n())
{
    allocate();

    for (label i=0; i<mRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i,j) = Mb(i,j);
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::~Matrix()
{
    if (v_)
    {
        delete[] v_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::clear()
{
    if (v_)
    {
        delete[] v_;
        v_ = nullptr;
    }

    mRows_ = 0;
    nCols_ = 0;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::transfer(Matrix<Form, Type>& M)
{
    clear();

    mRows_ = M.mRows_;
    M.mRows_ = 0;

    nCols_ = M.nCols_;
    M.nCols_ = 0;

    v_ = M.v_;
    M.v_ = nullptr;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::setSize(const label m, const label n)
{
    mType newMatrix(m, n, Zero);

    label minM = min(m, mRows_);
    label minN = min(n, nCols_);

    for (label i=0; i<minM; i++)
    {
        for (label j=0; j<minN; j++)
        {
            newMatrix(i, j) = (*this)(i, j);
        }
    }

    transfer(newMatrix);
}


template<class Form, class Type>
Form Foam::Matrix<Form, Type>::T() const
{
    const Matrix<Form, Type>& A = *this;
    Form At(n(), m());

    for (label i=0; i<m(); i++)
    {
        for (label j=0; j<n(); j++)
        {
            At(j, i) = A(i, j);
        }
    }

    return At;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Matrix<Form, Type>& M)
{
    if (this == &M)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (mRows_ != M.mRows_ || nCols_ != M.nCols_)
    {
        clear();
        mRows_ = M.mRows_;
        nCols_ = M.nCols_;
        allocate();
    }

    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = M.v_[i];
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
void Foam::Matrix<Form, Type>::operator=
(
    const ConstMatrixBlock<MatrixType>& Mb
)
{
    for (label i=0; i<mRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i,j) = Mb(i,j);
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
void Foam::Matrix<Form, Type>::operator=
(
    const MatrixBlock<MatrixType>& Mb
)
{
    for (label i=0; i<mRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i,j) = Mb(i,j);
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Type& s)
{
    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = s;
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const zero)
{
    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
const Type& Foam::max(const Matrix<Form, Type>& M)
{
    const label mn = M.size();

    if (mn)
    {
        label curMaxI = 0;
        const Type* Mv = M.v();

        for (label i=1; i<mn; i++)
        {
            if (Mv[i] > Mv[curMaxI])
            {
                curMaxI = i;
            }
        }

        return Mv[curMaxI];
    }
    else
    {
        FatalErrorInFunction
            << "Matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return M(0, 0);
    }
}


template<class Form, class Type>
const Type& Foam::min(const Matrix<Form, Type>& M)
{
    const label mn = M.size();

    if (mn)
    {
        label curMinI = 0;
        const Type* Mv = M.v();

        for (label i=1; i<mn; i++)
        {
            if (Mv[i] < Mv[curMinI])
            {
                curMinI = i;
            }
        }

        return Mv[curMinI];
    }
    else
    {
        FatalErrorInFunction
            << "Matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return M(0, 0);
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& M)
{
    Form nM(M.m(), M.n());

    if (M.m() && M.n())
    {
        Type* nMv = nM.v();
        const Type* Mv = M.v();

        const label mn = M.size();
        for (label i=0; i<mn; i++)
        {
            nMv[i] = -Mv[i];
        }
    }

    return nM;
}


template<class Form, class Type>
Form Foam::operator+(const Matrix<Form, Type>& A, const Matrix<Form, Type>& B)
{
    if (A.m() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different numbers of rows: "
            << A.m() << ", " << B.m()
            << abort(FatalError);
    }

    if (A.n() != B.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different numbers of columns: "
            << A.n() << ", " << B.n()
            << abort(FatalError);
    }

    Form AB(A.m(), A.n());

    Type* ABv = AB.v();
    const Type* Av = A.v();
    const Type* Bv = B.v();

    const label mn = A.size();
    for (label i=0; i<mn; i++)
    {
        ABv[i] = Av[i] + Bv[i];
    }

    return AB;
}


template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& A, const Matrix<Form, Type>& B)
{
    if (A.m() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different numbers of rows: "
            << A.m() << ", " << B.m()
            << abort(FatalError);
    }

    if (A.n() != B.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different numbers of columns: "
            << A.n() << ", " << B.n()
            << abort(FatalError);
    }

    Form AB(A.m(), A.n());

    Type* ABv = AB.v();
    const Type* Av = A.v();
    const Type* Bv = B.v();

    const label mn = A.size();
    for (label i=0; i<mn; i++)
    {
        ABv[i] = Av[i] - Bv[i];
    }

    return AB;
}


template<class Form, class Type>
Form Foam::operator*(const scalar s, const Matrix<Form, Type>& M)
{
    Form sM(M.m(), M.n());

    if (M.m() && M.n())
    {
        Type* sMv = sM.v();
        const Type* Mv = M.v();

        const label mn = M.size();
        for (label i=0; i<mn; i++)
        {
            sMv[i] = s*Mv[i];
        }
    }

    return sM;
}


template<class Form, class Type>
Form Foam::operator*(const Matrix<Form, Type>& M, const scalar s)
{
    Form sM(M.m(), M.n());

    if (M.m() && M.n())
    {
        Type* sMv = sM.v();
        const Type* Mv = M.v();

        const label mn = M.size();
        for (label i=0; i<mn; i++)
        {
            sMv[i] = Mv[i]*s;
        }
    }

    return sM;
}


template<class Form, class Type>
Form Foam::operator/(const Matrix<Form, Type>& M, const scalar s)
{
    Form sM(M.m(), M.n());

    if (M.m() && M.n())
    {
        Type* sMv = sM.v();
        const Type* Mv = M.v();

        const label mn = M.size();
        for (label i=0; i<mn; i++)
        {
            sMv[i] = Mv[i]/s;
        }
    }

    return sM;
}


template<class Form1, class Form2, class Type>
typename Foam::typeOfInnerProduct<Type, Form1, Form2>::type
Foam::operator*
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B
)
{
    if (A.n() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrices:" << nl
            << "Matrix A : " << A.m() << " x " << A.n() << nl
            << "Matrix B : " << B.m() << " x " << B.n() << nl
            << "In order to multiply matrices, columns of A must equal "
            << "rows of B"
            << abort(FatalError);
    }

    typename typeOfInnerProduct<Type, Form1, Form2>::type AB
    (
        A.m(),
        B.n(),
        Zero
    );

    for (label i=0; i<AB.m(); i++)
    {
        for (label j=0; j<AB.n(); j++)
        {
            for (label k=0; k<B.m(); k++)
            {
                AB(i, j) += A(i, k)*B(k, j);
            }
        }
    }

    return AB;
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::operator*
(
    const Matrix<Form, Type>& M,
    const Field<Type>& f
)
{
    if (M.n() != f.size())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrix and field:" << nl
            << "Matrix : " << M.m() << " x " << M.n() << nl
            << "Field : " << f.size() << " rows" << nl
            << "In order to multiply a matrix M and field f, "
               "columns of M must equal rows of f"
            << abort(FatalError);
    }

    tmp<Field<Type>> tMf(new Field<Type>(M.m(), Zero));
    Field<Type>& Mf = tMf.ref();

    for (label i=0; i<M.m(); i++)
    {
        for (label j=0; j<M.n(); j++)
        {
            Mf[i] += M(i, j)*f[j];
        }
    }

    return tMf;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "MatrixIO.C"

// ************************************************************************* //
