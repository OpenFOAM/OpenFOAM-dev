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

#include "LLTMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LLTMatrix<Type>::LLTMatrix()
{}


template<class Type>
Foam::LLTMatrix<Type>::LLTMatrix(const SquareMatrix<Type>& M)
{
    decompose(M);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LLTMatrix<Type>::decompose(const SquareMatrix<Type>& M)
{
    SquareMatrix<Type>& LLT = *this;

    // Initialize the LLT decomposition matrix to M
    LLT = M;

    const label m = LLT.m();

    for (label i=0; i<m; i++)
    {
        for (label j=0; j<m; j++)
        {
            if (j > i)
            {
                LLT(i, j) = Zero;
                continue;
            }

            Type sum = LLT(i, j);

            for (label k=0; k<j; k++)
            {
                sum -= LLT(i, k)*LLT(j, k);
            }

            if (i > j)
            {
                LLT(i, j) = sum/LLT(j, j);
            }
            else if (sum > 0)
            {
                LLT(i, i) = sqrt(sum);
            }
            else
            {
                FatalErrorInFunction
                    << "Cholesky decomposition failed, "
                       "matrix is not symmetric positive definite"
                    << abort(FatalError);
            }
        }
    }
}


template<class Type>
void Foam::LLTMatrix<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& source
) const
{
    // If x and source are different initialize x = source
    if (&x != &source)
    {
        x = source;
    }

    const SquareMatrix<Type>& LLT = *this;
    const label m = LLT.m();

    for (label i=0; i<m; i++)
    {
        Type sum = source[i];

        for (label j=0; j<i; j++)
        {
            sum = sum - LLT(i, j)*x[j];
        }

        x[i] = sum/LLT(i, i);
    }

    for (int i=m - 1; i >= 0; i--)
    {
        Type sum = x[i];

        for (label j=i + 1; j<m; j++)
        {
            sum = sum - LLT(j, i)*x[j];
        }

        x[i] = sum/LLT(i, i);
    }

}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LLTMatrix<Type>::solve
(
    const Field<Type>& source
) const
{
    tmp<Field<Type>> tx(new Field<Type>(this->m()));
    Field<Type>& x = tx.ref();

    solve(x, source);

    return tx;
}


// ************************************************************************* //
