/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix()
:
    Matrix<SymmetricSquareMatrix<Type>, Type>()
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix(const label n)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n)
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label m,
    const label n
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(m, n)
{
    if (m != n)
    {
        FatalErrorIn
        (
            "SymmetricSquareMatrix<Type>::SymmetricSquareMatrix"
            "(const label m, const label n)"
        )   << "m != n for constructing a symmetric square matrix"
            << exit(FatalError);
    }
}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label m,
    const label n,
    const Type& t
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(m, n, t)
{
    if (m != n)
    {
        FatalErrorIn
        (
            "SymmetricSquareMatrix<Type>::SymmetricSquareMatrix"
            "(const label m, const label n, const Type&)"
        )   << "m != n for constructing a symmetric square matrix"
            << exit(FatalError);
    }
}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix(Istream& is)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(is)
{}


template<class Type>
inline Foam::autoPtr<Foam::SymmetricSquareMatrix<Type> >
Foam::SymmetricSquareMatrix<Type>::clone() const
{
    return autoPtr<SymmetricSquareMatrix<Type> >
    (
        new SymmetricSquareMatrix<Type>(*this)
    );
}


template<class Type>
inline Type& Foam::SymmetricSquareMatrix<Type>::operator()
(
    const label r,
    const label c
)
{
    if (r > c)
    {
        return this->operator[](r)[c];
    }
    else
    {
        return this->operator[](c)[r];
    }
}


template<class Type>
inline const Type& Foam::SymmetricSquareMatrix<Type>::operator()
(
    const label r,
    const label c
) const
{
    if (r > c)
    {
        return this->operator[](r)[c];
    }
    else
    {
        return this->operator[](c)[r];
    }
}


// ************************************************************************* //
