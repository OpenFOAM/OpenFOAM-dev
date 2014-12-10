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

#include "DILUSmoother.H"
#include "DILUPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DILUSmoother, 0);

    lduMatrix::smoother::addasymMatrixConstructorToTable<DILUSmoother>
        addDILUSmootherAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DILUSmoother::DILUSmoother
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
    rD_(matrix_.diag())
{
    DILUPreconditioner::calcReciprocalD(rD_, matrix_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DILUSmoother::smooth
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    const scalar* const __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix_.lower().begin();

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

        register label nFaces = matrix_.upper().size();
        for (register label face=0; face<nFaces; face++)
        {
            register label u = uPtr[face];
            rAPtr[u] -= rDPtr[u]*lowerPtr[face]*rAPtr[lPtr[face]];
        }

        register label nFacesM1 = nFaces - 1;
        for (register label face=nFacesM1; face>=0; face--)
        {
            register label l = lPtr[face];
            rAPtr[l] -= rDPtr[l]*upperPtr[face]*rAPtr[uPtr[face]];
        }

        psi += rA;
    }
}


// ************************************************************************* //
