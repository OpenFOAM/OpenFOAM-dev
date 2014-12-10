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

Description
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::Amul
(
    scalarField& Apsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    scalar* __restrict__ ApsiPtr = Apsi.begin();

    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    register const label nCells = diag().size();
    for (register label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }


    register const label nFaces = upper().size();

    for (register label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    tpsi.clear();
}


void Foam::lduMatrix::Tmul
(
    scalarField& Tpsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    scalar* __restrict__ TpsiPtr = Tpsi.begin();

    const scalarField& psi = tpsi();
    const scalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ lowerPtr = lower().begin();
    const scalar* const __restrict__ upperPtr = upper().begin();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    register const label nCells = diag().size();
    for (register label cell=0; cell<nCells; cell++)
    {
        TpsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    register const label nFaces = upper().size();
    for (register label face=0; face<nFaces; face++)
    {
        TpsiPtr[uPtr[face]] += upperPtr[face]*psiPtr[lPtr[face]];
        TpsiPtr[lPtr[face]] += lowerPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    tpsi.clear();
}


void Foam::lduMatrix::sumA
(
    scalarField& sumA,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    scalar* __restrict__ sumAPtr = sumA.begin();

    const scalar* __restrict__ diagPtr = diag().begin();

    const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* __restrict__ lowerPtr = lower().begin();
    const scalar* __restrict__ upperPtr = upper().begin();

    register const label nCells = diag().size();
    register const label nFaces = upper().size();

    for (register label cell=0; cell<nCells; cell++)
    {
        sumAPtr[cell] = diagPtr[cell];
    }

    for (register label face=0; face<nFaces; face++)
    {
        sumAPtr[uPtr[face]] += lowerPtr[face];
        sumAPtr[lPtr[face]] += upperPtr[face];
    }

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces, patchI)
    {
        if (interfaces.set(patchI))
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            const scalarField& pCoeffs = interfaceBouCoeffs[patchI];

            forAll(pa, face)
            {
                sumAPtr[pa[face]] -= pCoeffs[face];
            }
        }
    }
}


void Foam::lduMatrix::residual
(
    scalarField& rA,
    const scalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    scalar* __restrict__ rAPtr = rA.begin();

    const scalar* const __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = diag().begin();
    const scalar* const __restrict__ sourcePtr = source.begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    FieldField<Field, scalar> mBouCoeffs(interfaceBouCoeffs.size());

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces.set(patchi))
        {
            mBouCoeffs.set(patchi, -interfaceBouCoeffs[patchi]);
        }
    }

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        mBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );

    register const label nCells = diag().size();
    for (register label cell=0; cell<nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - diagPtr[cell]*psiPtr[cell];
    }


    register const label nFaces = upper().size();

    for (register label face=0; face<nFaces; face++)
    {
        rAPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
        rAPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        mBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );
}


Foam::tmp<Foam::scalarField> Foam::lduMatrix::residual
(
    const scalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    tmp<scalarField> trA(new scalarField(psi.size()));
    residual(trA(), psi, source, interfaceBouCoeffs, interfaces, cmpt);
    return trA;
}


Foam::tmp<Foam::scalarField > Foam::lduMatrix::H1() const
{
    tmp<scalarField > tH1
    (
        new scalarField(lduAddr().size(), 0.0)
    );

    if (lowerPtr_ || upperPtr_)
    {
        scalarField& H1_ = tH1();

        scalar* __restrict__ H1Ptr = H1_.begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* __restrict__ lowerPtr = lower().begin();
        const scalar* __restrict__ upperPtr = upper().begin();

        register const label nFaces = upper().size();

        for (register label face=0; face<nFaces; face++)
        {
            H1Ptr[uPtr[face]] -= lowerPtr[face];
            H1Ptr[lPtr[face]] -= upperPtr[face];
        }
    }

    return tH1;
}


// ************************************************************************* //
