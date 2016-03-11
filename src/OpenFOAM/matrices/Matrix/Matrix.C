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
        v_ = new Type*[nRows_];
        v_[0] = new Type[nRows_*nCols_];

        for (label i=1; i<nRows_; i++)
        {
            v_[i] = v_[i-1] + nCols_;
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::~Matrix()
{
    if (v_)
    {
        delete[] (v_[0]);
        delete[] v_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label n, const label m)
:
    nRows_(n),
    nCols_(m),
    v_(NULL)
{
    if (nRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "bad n, m " << nRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label n, const label m, const Type& a)
:
    nRows_(n),
    nCols_(m),
    v_(NULL)
{
    if (nRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "bad n, m " << nRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        Type* v = v_[0];

        label mn = nRows_*nCols_;

        for (label i=0; i<mn; i++)
        {
            v[i] = a;
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
        Type* v = v_[0];
        const Type* av = a.v_[0];

        label mn = nRows_*nCols_;
        for (label i=0; i<mn; i++)
        {
            v[i] = av[i];
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::clear()
{
    if (v_)
    {
        delete[] (v_[0]);
        delete[] v_;
    }
    nRows_ = 0;
    nCols_ = 0;
    v_ = NULL;
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
Form Foam::Matrix<Form, Type>::T() const
{
    const Matrix<Form, Type>& A = *this;
    Form At(n(), m());

    for (label i=0; i<m(); i++)
    {
        for (label j=0; j<n(); j++)
        {
            At[j][i] = A[i][j];
        }
    }

    return At;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Type& t)
{
    if (v_)
    {
        Type* v = v_[0];

        label mn = nRows_*nCols_;
        for (label i=0; i<mn; i++)
        {
            v[i] = t;
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Matrix<Form, Type>& a)
{
    if (this == &a)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
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
        Type* v = v_[0];
        const Type* av = a.v_[0];

        label mn = nRows_*nCols_;
        for (label i=0; i<mn; i++)
        {
            v[i] = av[i];
        }
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
const Type& Foam::max(const Matrix<Form, Type>& a)
{
    label mn = a.m()*a.n();

    if (mn)
    {
        label curMaxI = 0;
        const Type* v = a[0];

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
            << "matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return a[0][0];
    }
}


template<class Form, class Type>
const Type& Foam::min(const Matrix<Form, Type>& a)
{
    label mn = a.m()*a.n();

    if (mn)
    {
        label curMinI = 0;
        const Type* v = a[0];

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
            << "matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return a[0][0];
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& a)
{
    Form na(a.m(), a.n());

    if (a.m() && a.n())
    {
        Type* nav = na[0];
        const Type* av = a[0];

        label mn = a.m()*a.n();
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
            << "attempted add matrices with different number of rows: "
            << a.m() << ", " << b.m()
            << abort(FatalError);
    }

    if (a.n() != b.n())
    {
        FatalErrorInFunction
            << "attempted add matrices with different number of columns: "
            << a.n() << ", " << b.n()
            << abort(FatalError);
    }

    Form ab(a.m(), a.n());

    Type* abv = ab[0];
    const Type* av = a[0];
    const Type* bv = b[0];

    label mn = a.m()*a.n();
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
            << "attempted add matrices with different number of rows: "
            << a.m() << ", " << b.m()
            << abort(FatalError);
    }

    if (a.n() != b.n())
    {
        FatalErrorInFunction
            << "attempted add matrices with different number of columns: "
            << a.n() << ", " << b.n()
            << abort(FatalError);
    }

    Form ab(a.m(), a.n());

    Type* abv = ab[0];
    const Type* av = a[0];
    const Type* bv = b[0];

    label mn = a.m()*a.n();
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
        Type* sav = sa[0];
        const Type* av = a[0];

        label mn = a.m()*a.n();
        for (label i=0; i<mn; i++)
        {
            sav[i] = s*av[i];
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
            << "attempted to multiply incompatible matrices:" << nl
            << "Matrix A : " << a.m() << " rows, " << a.n() << " columns" << nl
            << "Matrix B : " << b.m() << " rows, " << b.n() << " columns" << nl
            << "In order to multiply matrices, columns of A must equal "
            << "rows of B"
            << abort(FatalError);
    }

    Form ab(a.m(), b.n(), scalar(0));

    for (label i = 0; i < ab.m(); i++)
    {
        for (label j = 0; j < ab.n(); j++)
        {
            for (label l = 0; l < b.m(); l++)
            {
                ab[i][j] += a[i][l]*b[l][j];
            }
        }
    }

    return ab;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "MatrixIO.C"

// ************************************************************************* //
