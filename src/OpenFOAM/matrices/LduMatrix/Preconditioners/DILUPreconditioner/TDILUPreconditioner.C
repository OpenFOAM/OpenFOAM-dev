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

#include "TDILUPreconditioner.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::TDILUPreconditioner<Type, DType, LUType>::TDILUPreconditioner
(
    const typename LduMatrix<Type, DType, LUType>::solver& sol,
    const dictionary&
)
:
    LduMatrix<Type, DType, LUType>::preconditioner(sol),
    rD_(sol.matrix().diag())
{
    calcInvD(rD_, sol.matrix());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::TDILUPreconditioner<Type, DType, LUType>::calcInvD
(
    Field<DType>& rD,
    const LduMatrix<Type, DType, LUType>& matrix
)
{
    DType* __restrict__ rDPtr = rD.begin();

    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();

    const LUType* const __restrict__ upperPtr = matrix.upper().begin();
    const LUType* const __restrict__ lowerPtr = matrix.lower().begin();

    label nFaces = matrix.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -=
            dot(dot(upperPtr[face], lowerPtr[face]), inv(rDPtr[lPtr[face]]));
    }


    // Calculate the reciprocal of the preconditioned diagonal
    label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = inv(rDPtr[cell]);
    }
}


template<class Type, class DType, class LUType>
void Foam::TDILUPreconditioner<Type, DType, LUType>::precondition
(
    Field<Type>& wA,
    const Field<Type>& rA
) const
{
    Type* __restrict__ wAPtr = wA.begin();
    const Type* __restrict__ rAPtr = rA.begin();
    const DType* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        this->solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        this->solver_.matrix().lduAddr().lowerAddr().begin();
    const label* const __restrict__ losortPtr =
        this->solver_.matrix().lduAddr().losortAddr().begin();

    const LUType* const __restrict__ upperPtr =
        this->solver_.matrix().upper().begin();
    const LUType* const __restrict__ lowerPtr =
        this->solver_.matrix().lower().begin();

    label nCells = wA.size();
    label nFaces = this->solver_.matrix().upper().size();
    label nFacesM1 = nFaces - 1;

    for (label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = dot(rDPtr[cell], rAPtr[cell]);
    }


    label sface;

    for (label face=0; face<nFaces; face++)
    {
        sface = losortPtr[face];
        wAPtr[uPtr[sface]] -=
            dot(rDPtr[uPtr[sface]], dot(lowerPtr[sface], wAPtr[lPtr[sface]]));
    }

    for (label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -=
            dot(rDPtr[lPtr[face]], dot(upperPtr[face], wAPtr[uPtr[face]]));
    }
}


template<class Type, class DType, class LUType>
void Foam::TDILUPreconditioner<Type, DType, LUType>::preconditionT
(
    Field<Type>& wT,
    const Field<Type>& rT
) const
{
    Type* __restrict__ wTPtr = wT.begin();
    const Type* __restrict__ rTPtr = rT.begin();
    const DType* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        this->solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        this->solver_.matrix().lduAddr().lowerAddr().begin();
    const label* const __restrict__ losortPtr =
        this->solver_.matrix().lduAddr().losortAddr().begin();

    const LUType* const __restrict__ upperPtr =
        this->solver_.matrix().upper().begin();
    const LUType* const __restrict__ lowerPtr =
        this->solver_.matrix().lower().begin();

    label nCells = wT.size();
    label nFaces = this->solver_.matrix().upper().size();
    label nFacesM1 = nFaces - 1;

    for (label cell=0; cell<nCells; cell++)
    {
        wTPtr[cell] = dot(rDPtr[cell], rTPtr[cell]);
    }

    for (label face=0; face<nFaces; face++)
    {
        wTPtr[uPtr[face]] -=
            dot(rDPtr[uPtr[face]], dot(upperPtr[face], wTPtr[lPtr[face]]));
    }


    label sface;

    for (label face=nFacesM1; face>=0; face--)
    {
        sface = losortPtr[face];
        wTPtr[lPtr[sface]] -=
            dot(rDPtr[lPtr[sface]], dot(lowerPtr[sface], wTPtr[uPtr[sface]]));
    }
}


// ************************************************************************* //
