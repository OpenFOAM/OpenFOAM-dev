/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    const label n,
    const zero
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, Zero)
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    const Identity<Type>
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, Zero)
{
    for (label i=0; i<n; i++)
    {
        this->operator()(i, i) = pTraits<Type>::one;
    }
}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    const Type& t
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, t)
{}


template<class Type>
template<class InputIterator>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    InputIterator first,
    InputIterator last
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, first, last)
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    const label n,
    std::initializer_list<Type> lst
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(n, n, lst)
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    std::initializer_list<std::initializer_list<Type>> lstLst
)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(lstLst)
{
    if (this->m() != this->n())
    {
        FatalErrorInFunction
            << "Attempt to construct a symmetric square matrix "
               "from a rectangular matrix "
            << this->m() << " x " << this->n() << nl
            << abort(FatalError);
    }
}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix
(
    std::initializer_list<Type> lst
)
:
    SymmetricSquareMatrix<Type>(sqrt(scalar(lst.size())), lst)
{}


template<class Type>
inline Foam::SymmetricSquareMatrix<Type>::SymmetricSquareMatrix(Istream& is)
:
    Matrix<SymmetricSquareMatrix<Type>, Type>(is)
{}


template<class Type>
inline Foam::autoPtr<Foam::SymmetricSquareMatrix<Type>>
Foam::SymmetricSquareMatrix<Type>::clone() const
{
    return autoPtr<SymmetricSquareMatrix<Type>>
    (
        new SymmetricSquareMatrix<Type>(*this)
    );
}


// ************************************************************************* //
