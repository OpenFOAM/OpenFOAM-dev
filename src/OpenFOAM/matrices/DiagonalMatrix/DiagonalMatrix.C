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

#include "DiagonalMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::DiagonalMatrix<Type>::DiagonalMatrix()
:
    List<Type>()
{}


template<class Type>
template<class Form>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const Matrix<Form, Type>& a)
:
    List<Type>(min(a.n(), a.m()))
{
    forAll(*this, i)
    {
        this->operator[](i) = a[i][i];
    }
}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label size)
:
    List<Type>(size)
{}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label size, const Type& val)
:
    List<Type>(size, val)
{}


template<class Type>
Foam::DiagonalMatrix<Type>& Foam::DiagonalMatrix<Type>::invert()
{
    forAll(*this, i)
    {
        Type x = this->operator[](i);
        if (mag(x) < VSMALL)
        {
            this->operator[](i) = Type(0);
        }
        else
        {
            this->operator[](i) = Type(1)/x;
        }
    }

    return this;
}


template<class Type>
Foam::DiagonalMatrix<Type> Foam::inv(const DiagonalMatrix<Type>& A)
{
    DiagonalMatrix<Type> Ainv = A;

    forAll(A, i)
    {
        Type x = A[i];
        if (mag(x) < VSMALL)
        {
            Ainv[i] = Type(0);
        }
        else
        {
            Ainv[i] = Type(1)/x;
        }
    }

    return Ainv;
}


// ************************************************************************* //
