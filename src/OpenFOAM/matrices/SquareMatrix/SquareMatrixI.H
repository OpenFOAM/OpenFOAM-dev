/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
inline Foam::SquareMatrix<Type>::SquareMatrix(const label m, const label n)
:
    Matrix<SquareMatrix<Type>, Type>(m, n)
{
    if (m != n)
    {
        FatalErrorIn
        (
            "SquareMatrix<Type>::SquareMatrix(const label m, const label n)"
        ) << "m != n for constructing a square matrix" << exit(FatalError);
    }
}

template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label m,
    const label n,
    const Type& t
)
:
    Matrix<SquareMatrix<Type>, Type>(m, n, t)
{
    if (m != n)
    {
        FatalErrorIn
        (
            "SquareMatrix<Type>::SquareMatrix"
            "(const label m, const label n, const Type&)"
        ) << "m != n for constructing a square matrix" << exit(FatalError);
    }
}

template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(Istream& is)
:
    Matrix<SquareMatrix<Type>, Type>(is)
{}

template<class Type>
inline Foam::autoPtr<Foam::SquareMatrix<Type> >
Foam::SquareMatrix<Type>::clone() const
{
    return autoPtr<SquareMatrix<Type> >(new SquareMatrix<Type>(*this));
}


// ************************************************************************* //
