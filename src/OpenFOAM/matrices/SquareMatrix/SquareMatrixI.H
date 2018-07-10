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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix()
:
    Matrix<SquareMatrix<Type>, Type>()
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(const label n)
:
    Matrix<SquareMatrix<Type>, Type>(n, n)
{}


template<class Type>
template<class MatrixType>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const ConstMatrixBlock<MatrixType>& block
)
:
    Matrix<SquareMatrix<Type>, Type>(block)
{}


template<class Type>
template<class MatrixType>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const MatrixBlock<MatrixType>& block
)
:
    Matrix<SquareMatrix<Type>, Type>(block)
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const zero
)
:
    Matrix<SquareMatrix<Type>, Type>(n, n, Zero)
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label m,
    const label n,
    const zero
)
:
    Matrix<SquareMatrix<Type>, Type>(m, n, Zero)
{
    if (m != n)
    {
        FatalErrorInFunction
            << "Attempt to construct a square matrix "
            << m << " x " << n << nl
            << abort(FatalError);
    }
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const Identity<Type>
)
:
    Matrix<SquareMatrix<Type>, Type>(n, n, Zero)
{
    for (label i=0; i<n; i++)
    {
        this->operator()(i, i) = Type(I);
    }
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const Type& t
)
:
    Matrix<SquareMatrix<Type>, Type>(n, n, t)
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const RectangularMatrix<Type>& RM
)
:
    Matrix<SquareMatrix<Type>, Type>(RM)
{
    if (this->m() != this->n())
    {
        FatalErrorInFunction
            << "Attempt to construct a square matrix from a rectangular matrix "
            << this->m() << " x " << this->n() << nl
            << abort(FatalError);
    }
}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(Istream& is)
:
    Matrix<SquareMatrix<Type>, Type>(is)
{}


template<class Type>
inline Foam::autoPtr<Foam::SquareMatrix<Type>>
Foam::SquareMatrix<Type>::clone() const
{
    return autoPtr<SquareMatrix<Type>>(new SquareMatrix<Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SquareMatrix<Type>::setSize(const label m)
{
    Matrix<SquareMatrix<Type>, Type>::setSize(m, m);
}


template<class Type>
inline void Foam::SquareMatrix<Type>::shallowResize(const label m)
{
    Matrix<SquareMatrix<Type>, Type>::shallowResize(m, m);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::SquareMatrix<Type>::operator=(const zero)
{
    Matrix<SquareMatrix<Type>, Type>::operator=(Zero);
}


template<class Type>
void Foam::SquareMatrix<Type>::operator=(const Identity<Type>)
{
    Matrix<SquareMatrix<Type>, Type>::operator=(Zero);
    for (label i=0; i<this->n(); i++)
    {
        this->operator()(i, i) = Type(I);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::SquareMatrix<Type> symmOuter
(
    const Field<Type>& f1,
    const Field<Type>& f2
)
{
    SquareMatrix<Type> f1f2T(f1.size());

    for (label i=0; i<f1f2T.m(); i++)
    {
        for (label j=0; j<f1f2T.n(); j++)
        {
            f1f2T(i, j) = f1[i]*f2[j];
        }
    }

    return f1f2T;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
