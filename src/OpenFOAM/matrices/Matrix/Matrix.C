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

\*---------------------------------------------------------------------------*/

#include "Matrix.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::allocate()
{
    if (nRows_ && nCols_)
    {
        v_ = new Type[size()];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n)
:
    nRows_(m),
    nCols_(n),
    v_(NULL)
{
    if (nRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Bad m, n " << nRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const zero)
:
    nRows_(m),
    nCols_(n),
    v_(NULL)
{
    if (nRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Bad m, n " << nRows_ << ", " << nCols_
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
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const Type& a)
:
    nRows_(m),
    nCols_(n),
    v_(NULL)
{
    if (nRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Bad m, n " << nRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = a;
        }
    }
}


template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix(const mType::Block& block)
:
    nRows_(block.m()),
    nCols_(block.n())
{
    allocate();

    for (label i=0; i<nRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i,j) = block(i,j);
        }
    }
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form, Type>& a)
:
    nRows_(a.nRows_),
    nCols_(a.nCols_),
    v_(NULL)
{
    if (a.v_)
    {
        allocate();

        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = a.v_[i];
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
        v_ = NULL;
    }

    nRows_ = 0;
    nCols_ = 0;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::transfer(Matrix<Form, Type>& a)
{
    clear();

    nRows_ = a.nRows_;
    a.nRows_ = 0;

    nCols_ = a.nCols_;
    a.nCols_ = 0;

    v_ = a.v_;
    a.v_ = NULL;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::setSize(const label m, const label n)
{
    mType newMatrix(m, n, Zero);

    label minM = min(m, nRows_);
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


template<class Form, class Type>
Foam::Matrix<Form, Type>::ConstBlock::operator Field<cmptType>() const
{
    if (nCols_ != 1)
    {
        FatalErrorInFunction
            << "Number of columns " << nCols_ << " != 1"
            << abort(FatalError);
    }

    Field<cmptType> f(nRows_);

    forAll(f, i)
    {
        f[i] = operator()(i, 0);
    }

    return f;
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Block::operator Field<cmptType>() const
{
    if (nCols_ != 1)
    {
        FatalErrorInFunction
            << "Number of columns " << nCols_ << " != 1"
            << abort(FatalError);
    }

    Field<cmptType> f(nRows_);

    forAll(f, i)
    {
        f[i] = operator()(i, 0);
    }

    return f;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Matrix<Form, Type>& a)
{
    if (this == &a)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (nRows_ != a.nRows_ || nCols_ != a.nCols_)
    {
        clear();
        nRows_ = a.nRows_;
        nCols_ = a.nCols_;
        allocate();
    }

    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = a.v_[i];
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const mType::Block& block)
{
    for (label i=0; i<nRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i,j) = block(i,j);
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::Block::operator=(const Block& block)
{
    if (this != &block)
    {
        if (nRows_ != block.m() || nCols_ != block.n())
        {
            FatalErrorInFunction
                << "Attempt to assign blocks of different sizes: "
                << nRows_ << "x" << nCols_ << " != "
                << block.m() << "x" << block.n()
                << abort(FatalError);
        }

        for (label i=0; i<nRows_; i++)
        {
            for (label j=0; j<nCols_; j++)
            {
                (*this)(i, j) = block(i, j);
            }
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::Block::operator=(const ConstBlock& block)
{
    if (this != &block)
    {
        if (nRows_ != block.m() || nCols_ != block.n())
        {
            FatalErrorInFunction
                << "Attempt to assign blocks of different sizes: "
                << nRows_ << "x" << nCols_ << " != "
                << block.m() << "x" << block.n()
                << abort(FatalError);
        }

        for (label i=0; i<nRows_; i++)
        {
            for (label j=0; j<nCols_; j++)
            {
                (*this)(i, j) = block(i, j);
            }
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::Block::operator=(const mType& block)
{
    if (nRows_ != block.m() || nCols_ != block.n())
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << nRows_ << "x" << nCols_ << " != "
            << block.m() << "x" << block.n()
            << abort(FatalError);
    }

    for (label i=0; i<nRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i, j) = block(i, j);
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Type& t)
{
    if (v_)
    {
        const label mn = size();
        for (label i=0; i<mn; i++)
        {
            v_[i] = t;
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


template<class Form, class Type>
template
<
    template<class, Foam::direction, Foam::direction> class MSBlock,
    class SubTensor,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
void Foam::Matrix<Form, Type>::Block::operator=
(
    const MSBlock<SubTensor, BRowStart, BColStart>& block
)
{
    if (nRows_ != block.nRows || nCols_ != block.nCols)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << nRows_ << "x" << nCols_ << " != "
            << block.nRows << "x" << block.nCols
            << abort(FatalError);
    }

    for (direction i=0; i<nRows_; ++i)
    {
        for (direction j=0; j<nCols_; ++j)
        {
            operator()(i, j) = block(i, j);
        }
    }
}


template<class Form, class Type>
template
<
    template<class, Foam::direction> class VSBlock,
    class SubVector,
    Foam::direction BStart
>
void Foam::Matrix<Form, Type>::Block::operator=
(
    const VSBlock<SubVector, BStart>& block
)
{
    if (nRows_ != block.nComponents || nCols_ != 1)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << nRows_ << "x" << nCols_ << " != "
            << block.nComponents << "x" << 1
            << abort(FatalError);
    }

    for (direction i=0; i<nRows_; ++i)
    {
        operator()(i, 0) = block[i];
    }
}


template<class Form, class Type>
template<class MSForm, Foam::direction Nrows, Foam::direction Ncols>
void Foam::Matrix<Form, Type>::Block::operator=
(
    const MatrixSpace<MSForm, cmptType, Nrows, Ncols>& ms
)
{
    if (nRows_ != Nrows || nCols_ != Ncols)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << nRows_ << "x" << nCols_ << " != "
            << Nrows << "x" << Ncols
            << abort(FatalError);
    }

    for (label i=0; i<nRows_; i++)
    {
        for (label j=0; j<nCols_; j++)
        {
            (*this)(i, j) = ms(i, j);
        }
    }
}


template<class Form, class Type>
template<class VSForm, Foam::direction Ncmpts>
void Foam::Matrix<Form, Type>::Block::operator=
(
    const VectorSpace<VSForm, cmptType, Ncmpts>& ms
)
{
    if (nRows_ != Ncmpts || nCols_ != 1)
    {
        FatalErrorInFunction
            << "Attempt to assign blocks of different sizes: "
            << nRows_ << "x" << nCols_ << " != "
            << Ncmpts << "x" << 1
            << abort(FatalError);
    }

    for (direction i=0; i<Ncmpts; ++i)
    {
        operator()(i, 0) = ms[i];
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::Block::operator=(const Field<cmptType>& f)
{
    if (nRows_ != f.size() || nCols_ != 1)
    {
        FatalErrorInFunction
            << "Error: cannot assign blocks of different size (left is "
            << nRows_ << "x" << nCols_ << " != "
            << f.size() << "x" << 1
            << abort(FatalError);
    }

    forAll(f, i)
    {
        operator()(i, 0) = f[i];
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
const Type& Foam::max(const Matrix<Form, Type>& a)
{
    const label mn = a.size();

    if (mn)
    {
        label curMaxI = 0;
        const Type* v = a.v();

        for (label i=1; i<mn; i++)
        {
            if (v[i] > v[curMaxI])
            {
                curMaxI = i;
            }
        }

        return v[curMaxI];
    }
    else
    {
        FatalErrorInFunction
            << "Matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return a(0, 0);
    }
}


template<class Form, class Type>
const Type& Foam::min(const Matrix<Form, Type>& a)
{
    const label mn = a.size();

    if (mn)
    {
        label curMinI = 0;
        const Type* v = a.v();

        for (label i=1; i<mn; i++)
        {
            if (v[i] < v[curMinI])
            {
                curMinI = i;
            }
        }

        return v[curMinI];
    }
    else
    {
        FatalErrorInFunction
            << "Matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return a(0, 0);
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& a)
{
    Form na(a.m(), a.n());

    if (a.m() && a.n())
    {
        Type* nav = na.v();
        const Type* av = a.v();

        const label mn = a.size();
        for (label i=0; i<mn; i++)
        {
            nav[i] = -av[i];
        }
    }

    return na;
}


template<class Form, class Type>
Form Foam::operator+(const Matrix<Form, Type>& a, const Matrix<Form, Type>& b)
{
    if (a.m() != b.m())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different number of rows: "
            << a.m() << ", " << b.m()
            << abort(FatalError);
    }

    if (a.n() != b.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different number of columns: "
            << a.n() << ", " << b.n()
            << abort(FatalError);
    }

    Form ab(a.m(), a.n());

    Type* abv = ab.v();
    const Type* av = a.v();
    const Type* bv = b.v();

    const label mn = a.size();
    for (label i=0; i<mn; i++)
    {
        abv[i] = av[i] + bv[i];
    }

    return ab;
}


template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& a, const Matrix<Form, Type>& b)
{
    if (a.m() != b.m())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different number of rows: "
            << a.m() << ", " << b.m()
            << abort(FatalError);
    }

    if (a.n() != b.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different number of columns: "
            << a.n() << ", " << b.n()
            << abort(FatalError);
    }

    Form ab(a.m(), a.n());

    Type* abv = ab.v();
    const Type* av = a.v();
    const Type* bv = b.v();

    const label mn = a.size();
    for (label i=0; i<mn; i++)
    {
        abv[i] = av[i] - bv[i];
    }

    return ab;
}


template<class Form, class Type>
Form Foam::operator*(const scalar s, const Matrix<Form, Type>& a)
{
    Form sa(a.m(), a.n());

    if (a.m() && a.n())
    {
        Type* sav = sa.v();
        const Type* av = a.v();

        const label mn = a.size();
        for (label i=0; i<mn; i++)
        {
            sav[i] = s*av[i];
        }
    }

    return sa;
}


template<class Form, class Type>
Form Foam::operator*(const Matrix<Form, Type>& a, const scalar s)
{
    Form sa(a.m(), a.n());

    if (a.m() && a.n())
    {
        Type* sav = sa.v();
        const Type* av = a.v();

        const label mn = a.size();
        for (label i=0; i<mn; i++)
        {
            sav[i] = av[i]*s;
        }
    }

    return sa;
}


template<class Form, class Type>
Form Foam::operator/(const Matrix<Form, Type>& a, const scalar s)
{
    Form sa(a.m(), a.n());

    if (a.m() && a.n())
    {
        Type* sav = sa.v();
        const Type* av = a.v();

        const label mn = a.size();
        for (label i=0; i<mn; i++)
        {
            sav[i] = av[i]/s;
        }
    }

    return sa;
}


template<class Form, class Type>
Form Foam::operator*(const Matrix<Form, Type>& a, const Matrix<Form, Type>& b)
{
    if (a.n() != b.m())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrices:" << nl
            << "Matrix A : " << a.m() << " rows, " << a.n() << " columns" << nl
            << "Matrix B : " << b.m() << " rows, " << b.n() << " columns" << nl
            << "In order to multiply matrices, columns of A must equal "
            << "rows of B"
            << abort(FatalError);
    }

    Form ab(a.m(), b.n(), scalar(0));

    for (label i=0; i<ab.m(); i++)
    {
        for (label j=0; j<ab.n(); j++)
        {
            for (label k=0; k<b.m(); k++)
            {
                ab(i, j) += a(i, k)*b(k, j);
            }
        }
    }

    return ab;
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
            << "Matrix : " << M.m() << " rows, " << M.n() << " columns" << nl
            << "Field : " << f.size() << " rows" << nl
            << "In order to multiply a matrix M and field f, "
               "columns of M must equal rows of f"
            << abort(FatalError);
    }

    tmp<Field<Type>> tMf(new Field<Type>(f.size(), Zero));
    Field<Type>& Mf = tMf.ref();

    for (label i=0; i<M.m(); ++i)
    {
        for (label j=0; j<M.n(); ++j)
        {
            Mf[i] += M(i, j)*f[j];
        }
    }

    return tMf;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "MatrixIO.C"

// ************************************************************************* //
