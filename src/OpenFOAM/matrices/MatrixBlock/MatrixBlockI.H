/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

template<class MatrixType>
Foam::ConstMatrixBlock<MatrixType>::ConstMatrixBlock
(
    const MatrixType& matrix,
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
:
    matrix_(matrix),
    mRows_(m),
    nCols_(n),
    rowStart_(mStart),
    colStart_(nStart)
{
    #ifdef FULLDEBUG
    if
    (
        rowStart_ + mRows_ > matrix.m()
     || colStart_ + nCols_ > matrix.n()
    )
    {
        FatalErrorInFunction
            << "Block addresses outside matrix"
            << abort(FatalError);
    }
    #endif
}


template<class MatrixType>
Foam::MatrixBlock<MatrixType>::MatrixBlock
(
    MatrixType& matrix,
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
:
    matrix_(matrix),
    mRows_(m),
    nCols_(n),
    rowStart_(mStart),
    colStart_(nStart)
{
    #ifdef FULLDEBUG
    if
    (
        rowStart_ + mRows_ > matrix.m()
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

template<class MatrixType>
inline Foam::label Foam::ConstMatrixBlock<MatrixType>::m() const
{
    return mRows_;
}


template<class MatrixType>
inline Foam::label Foam::ConstMatrixBlock<MatrixType>::n() const
{
    return nCols_;
}


template<class MatrixType>
inline Foam::label Foam::MatrixBlock<MatrixType>::m() const
{
    return mRows_;
}


template<class MatrixType>
inline Foam::label Foam::MatrixBlock<MatrixType>::n() const
{
    return nCols_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class MatrixType>
inline const typename MatrixType::cmptType&
Foam::ConstMatrixBlock<MatrixType>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    if (i<0 || i>=mRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << mRows_-1
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


template<class MatrixType>
inline const typename MatrixType::cmptType&
Foam::MatrixBlock<MatrixType>::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    if (i<0 || i>=mRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << mRows_-1
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


template<class MatrixType>
inline typename MatrixType::cmptType&
Foam::MatrixBlock<MatrixType>::operator()
(
    const label i,
    const label j
)
{
    #ifdef FULLDEBUG
    if (i<0 || i>=mRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << mRows_-1
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


// ************************************************************************* //
