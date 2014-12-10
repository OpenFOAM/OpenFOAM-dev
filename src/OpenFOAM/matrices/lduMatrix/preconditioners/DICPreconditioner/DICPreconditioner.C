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

#include "DICPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DICPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<DICPreconditioner>
        addDICPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DICPreconditioner::DICPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag())
{
    calcReciprocalD(rD_, sol.matrix());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DICPreconditioner::calcReciprocalD
(
    scalarField& rD,
    const lduMatrix& matrix
)
{
    scalar* __restrict__ rDPtr = rD.begin();

    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    // Calculate the DIC diagonal
    register const label nFaces = matrix.upper().size();
    for (register label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rDPtr[lPtr[face]];
    }


    // Calculate the reciprocal of the preconditioned diagonal
    register const label nCells = rD.size();

    for (register label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }
}


void Foam::DICPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction
) const
{
    scalar* __restrict__ wAPtr = wA.begin();
    const scalar* __restrict__ rAPtr = rA.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().begin();

    register label nCells = wA.size();
    register label nFaces = solver_.matrix().upper().size();
    register label nFacesM1 = nFaces - 1;

    for (register label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    for (register label face=0; face<nFaces; face++)
    {
        wAPtr[uPtr[face]] -= rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
    }

    for (register label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -= rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
    }
}


// ************************************************************************* //
