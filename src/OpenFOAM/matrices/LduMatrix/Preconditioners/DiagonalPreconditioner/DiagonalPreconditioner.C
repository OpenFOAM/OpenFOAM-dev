/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "DiagonalPreconditioner.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::DiagonalPreconditioner<Type, DType, LUType>::DiagonalPreconditioner
(
    const typename LduMatrix<Type, DType, LUType>::solver& sol,
    const dictionary&
)
:
    LduMatrix<Type, DType, LUType>::preconditioner(sol),
    rD(sol.matrix().diag().size())
{
    DType* __restrict__ rDPtr = rD.begin();
    const DType* __restrict__ DPtr = this->solver_.matrix().diag().begin();

    register label nCells = rD.size();

    // Generate inverse (reciprocal for scalar) diagonal
    for (register label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = inv(DPtr[cell]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::DiagonalPreconditioner<Type, DType, LUType>::read(const dictionary&)
{}


template<class Type, class DType, class LUType>
void Foam::DiagonalPreconditioner<Type, DType, LUType>::precondition
(
    Field<Type>& wA,
    const Field<Type>& rA
) const
{
    Type* __restrict__ wAPtr = wA.begin();
    const Type* __restrict__ rAPtr = rA.begin();
    const DType* __restrict__ rDPtr = rD.begin();

    register label nCells = wA.size();

    for (register label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = dot(rDPtr[cell], rAPtr[cell]);
    }
}


// ************************************************************************* //
