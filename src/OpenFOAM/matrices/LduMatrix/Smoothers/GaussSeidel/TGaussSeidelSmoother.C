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

#include "TGaussSeidelSmoother.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::TGaussSeidelSmoother<Type, DType, LUType>::TGaussSeidelSmoother
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix
)
:
    LduMatrix<Type, DType, LUType>::smoother
    (
        fieldName,
        matrix
    ),
    rD_(matrix.diag().size())
{
    register const label nCells = matrix.diag().size();
    register const DType* const __restrict__ diagPtr = matrix.diag().begin();
    register DType* __restrict__ rDPtr = rD_.begin();

    for (register label cellI=0; cellI<nCells; cellI++)
    {
        rDPtr[cellI] = inv(diagPtr[cellI]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::TGaussSeidelSmoother<Type, DType, LUType>::smooth
(
    const word& fieldName_,
    Field<Type>& psi,
    const LduMatrix<Type, DType, LUType>& matrix_,
    const Field<DType>& rD_,
    const label nSweeps
)
{
    register Type* __restrict__ psiPtr = psi.begin();

    register const label nCells = psi.size();

    Field<Type> bPrime(nCells);
    register Type* __restrict__ bPrimePtr = bPrime.begin();

    register const DType* const __restrict__ rDPtr = rD_.begin();

    register const LUType* const __restrict__ upperPtr =
        matrix_.upper().begin();

    register const LUType* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    register const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    register const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update to add the contibution to the r.h.s.

    FieldField<Field, LUType> mBouCoeffs(matrix_.interfacesUpper().size());

    forAll(mBouCoeffs, patchi)
    {
        if (matrix_.interfaces().set(patchi))
        {
            mBouCoeffs.set(patchi, -matrix_.interfacesUpper()[patchi]);
        }
    }

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        bPrime = matrix_.source();

        matrix_.initMatrixInterfaces
        (
            mBouCoeffs,
            psi,
            bPrime
        );

        matrix_.updateMatrixInterfaces
        (
            mBouCoeffs,
            psi,
            bPrime
        );

        Type curPsi;
        register label fStart;
        register label fEnd = ownStartPtr[0];

        for (register label cellI=0; cellI<nCells; cellI++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = ownStartPtr[cellI + 1];

            // Get the accumulated neighbour side
            curPsi = bPrimePtr[cellI];

            // Accumulate the owner product side
            for (register label curFace=fStart; curFace<fEnd; curFace++)
            {
                curPsi -= dot(upperPtr[curFace], psiPtr[uPtr[curFace]]);
            }

            // Finish current psi
            curPsi = dot(rDPtr[cellI], curPsi);

            // Distribute the neighbour side using current psi
            for (register label curFace=fStart; curFace<fEnd; curFace++)
            {
                bPrimePtr[uPtr[curFace]] -= dot(lowerPtr[curFace], curPsi);
            }

            psiPtr[cellI] = curPsi;
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::TGaussSeidelSmoother<Type, DType, LUType>::smooth
(
    Field<Type>& psi,
    const label nSweeps
) const
{
    smooth(this->fieldName_, psi, this->matrix_, rD_, nSweeps);
}


// ************************************************************************* //
