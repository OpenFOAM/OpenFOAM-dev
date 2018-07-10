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

#include "FDICSmoother.H"
#include "FDICPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FDICSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<FDICSmoother>
        addFDICSmootherSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FDICSmoother::FDICSmoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    rD_(matrix_.diag()),
    rDuUpper_(matrix_.upper().size()),
    rDlUpper_(matrix_.upper().size())
{
    scalar* __restrict__ rDPtr = rD_.begin();
    scalar* __restrict__ rDuUpperPtr = rDuUpper_.begin();
    scalar* __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();

    label nCells = rD_.size();
    label nFaces = matrix_.upper().size();

    for (label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= sqr(upperPtr[face])/rDPtr[lPtr[face]];
    }

    // Generate reciprocal FDIC
    for (label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        rDuUpperPtr[face] = rDPtr[uPtr[face]]*upperPtr[face];
        rDlUpperPtr[face] = rDPtr[lPtr[face]]*upperPtr[face];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FDICSmoother::smooth
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    const scalar* const __restrict__ rDuUpperPtr = rDuUpper_.begin();
    const scalar* const __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();

    // Temporary storage for the residual
    scalarField rA(rD_.size());
    scalar* __restrict__ rAPtr = rA.begin();

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        matrix_.residual
        (
            rA,
            psi,
            source,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt
        );

        rA *= rD_;

        label nFaces = matrix_.upper().size();
        for (label face=0; face<nFaces; face++)
        {
            rAPtr[uPtr[face]] -= rDuUpperPtr[face]*rAPtr[lPtr[face]];
        }

        label nFacesM1 = nFaces - 1;
        for (label face=nFacesM1; face>=0; face--)
        {
            rAPtr[lPtr[face]] -= rDlUpperPtr[face]*rAPtr[uPtr[face]];
        }

        psi += rA;
    }
}


// ************************************************************************* //
