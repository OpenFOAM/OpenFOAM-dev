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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix()
:
    nRows_(0),
    nCols_(0),
    v_(NULL)
{}


template<class Form, class Type>
inline Foam::autoPtr<Foam::Matrix<Form, Type>> Foam::Matrix<Form, Type>::
clone() const
{
    return autoPtr<Matrix<Form, Type>>(new Matrix<Form, Type>(*this));
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::ConstBlock::ConstBlock
(
    const mType& matrix,
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
:
    matrix_(matrix),
    nRows_(m),
    nCols_(n),
    rowStart_(mStart),
    colStart_(nStart)
{
    #ifdef FULLDEBUG
    if
    (
        rowStart_ + nRows_ > matrix.m()
     || colStart_ + nCols_ > matrix.n()
    )
    {
        FatalErrorInFunction
            << "Block addresses outside matrix"
            << abort(FatalError);
    }
    #endif
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Block::Block
(
    mType& matrix,
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
:
    matrix_(matrix),
    nRows_(m),
    nCols_(n),
    rowStart_(mStart),
    colStart_(nStart)
{
    #ifdef FULLDEBUG
    if
    (
        rowStart_ + nRows_ > matrix.m()
     || colStart_ + nCols_ > matrix.n()
    )
    {
        FatalErrorInFunction
            << "Block addresses outside matrix"
            << abort(FatalError);
    }
    #endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
inline const Foam::Matrix<Form, Type>& Foam::Matrix<Form, Type>::null()
{
    return NullObjectRef<Matrix<Form, Type>>();
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::m() const
{
    return nRows_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::n() const
{
    return nCols_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::size() const
{
    return nRows_*nCols_;
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checki(const label i) const
{
    #ifdef FULLDEBUG
    if (!nRows_ || !nCols_)
    {
        FatalErrorInFunction
            << "Attempt to access element from empty matrix"
            << abort(FatalError);
    }
    else if (i<0 || i>=nRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << nRows_-1
            << abort(FatalError);
    }
    #endif
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checkj(const label j) const
{
    #ifdef FULLDEBUG
    if (!nRows_ || !nCols_)
    {
        FatalErrorInFunction
            << "Attempt to access element from empty matrix"
            << abort(FatalError);
    }
    else if (j<0 || j>=nCols_)
    {
        FatalErrorInFunction
            << "index " << j << " out of range 0 ... " << nCols_-1
            << abort(FatalError);
    }
    #endif
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::ConstBlock::m() const
{
    return nRows_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::ConstBlock::n() const
{
    return nCols_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::Block::m() const
{
    return nRows_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::Block::n() const
{
    return nCols_;
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::v() const
{
    return v_;
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::v()
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::ConstBlock
Foam::Matrix<Form, Type>::block
(
    const label m,
    const label n,
    const label mStart,
    const label nStart
) const
{
    return ConstBlock
    (
        *this,
        m,
        n,
        mStart,
        nStart
    );
}


template<class Form, class Type>
template<class VectorSpace>
inline typename Foam::Matrix<Form, Type>::ConstBlock
Foam::Matrix<Form, Type>::block
(
    const label mStart,
    const label nStart
) const
{
    return ConstBlock
    (
        *this,
        VectorSpace::nRows,
        VectorSpace::nCols,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::ConstBlock
Foam::Matrix<Form, Type>::col
(
    const label m,
    const label mStart
) const
{
    return ConstBlock
    (
        *this,
        m,
        1,
        mStart,
        0
    );
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::ConstBlock
Foam::Matrix<Form, Type>::col
(
    const label m,
    const label mStart,
    const label nStart
) const
{
    return ConstBlock
    (
        *this,
        m,
        1,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::Block
Foam::Matrix<Form, Type>::block
(
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
{
    return Block
    (
        *this,
        m,
        n,
        mStart,
        nStart
    );
}


template<class Form, class Type>
template<class VectorSpace>
inline typename Foam::Matrix<Form, Type>::Block
Foam::Matrix<Form, Type>::block(const label mStart, const label nStart)
{
    return Block
    (
        *this,
        VectorSpace::nRows,
        VectorSpace::nCols,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::Block
Foam::Matrix<Form, Type>::col(const label m, const label mStart)
{
    return Block
    (
        *this,
        m,
        1,
        mStart,
        0
    );
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::Block
Foam::Matrix<Form, Type>::col
(
    const label m,
    const label mStart,
    const label nStart
)
{
    return Block
    (
        *this,
        m,
        1,
        mStart,
        nStart
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::ConstBlock::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    if (i<0 || i>=nRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << nRows_-1
            << abort(FatalError);
    }
    if (j<0 || j>=nCols_)
    {
        FatalErrorInFunction
            << "Index " << j << " out of range 0 ... " << nCols_-1
            << abort(FatalError);
    }
    #endif

    return matrix_(i + rowStart_, j + colStart_);
}


template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::Block::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    if (i<0 || i>=nRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << nRows_-1
            << abort(FatalError);
    }
    if (j<0 || j>=nCols_)
    {
        FatalErrorInFunction
            << "Index " << j << " out of range 0 ... " << nCols_-1
            << abort(FatalError);
    }
    #endif

    return matrix_(i + rowStart_, j + colStart_);
}


template<class Form, class Type>
inline Type& Foam::Matrix<Form, Type>::Block::operator()
(
    const label i,
    const label j
)
{
    #ifdef FULLDEBUG
    if (i<0 || i>=nRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << nRows_-1
            << abort(FatalError);
    }
    if (j<0 || j>=nCols_)
    {
        FatalErrorInFunction
            << "Index " << j << " out of range 0 ... " << nCols_-1
            << abort(FatalError);
    }
    #endif

    return matrix_(i + rowStart_, j + colStart_);
}


template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::operator()
(
    const label i,
    const label j
) const
{
    checki(i);
    checkj(j);
    return v_[i*nCols_ + j];
}


template<class Form, class Type>
inline Type& Foam::Matrix<Form, Type>::operator()
(
    const label i,
    const label j
)
{
    checki(i);
    checkj(j);
    return v_[i*nCols_ + j];
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::operator[](const label i) const
{
    checki(i);
    return v_ + i*nCols_;
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::operator[](const label i)
{
    checki(i);
    return v_ + i*nCols_;
}


// ************************************************************************* //
