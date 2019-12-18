/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    const label nCells = matrix.diag().size();
    const DType* const __restrict__ diagPtr = matrix.diag().begin();
    DType* __restrict__ rDPtr = rD_.begin();

    for (label celli=0; celli<nCells; celli++)
    {
        rDPtr[celli] = inv(diagPtr[celli]);
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
    Type* __restrict__ psiPtr = psi.begin();

    const label nCells = psi.size();

    Field<Type> bPrime(nCells);
    Type* __restrict__ bPrimePtr = bPrime.begin();

    const DType* const __restrict__ rDPtr = rD_.begin();

    const LUType* const __restrict__ upperPtr =
        matrix_.upper().begin();

    const LUType* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update to add the contribution to the r.h.s.

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
        label fStart;
        label fEnd = ownStartPtr[0];

        for (label celli=0; celli<nCells; celli++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = ownStartPtr[celli + 1];

            // Get the accumulated neighbour side
            curPsi = bPrimePtr[celli];

            // Accumulate the owner product side
            for (label curFace=fStart; curFace<fEnd; curFace++)
            {
                curPsi -= dot(upperPtr[curFace], psiPtr[uPtr[curFace]]);
            }

            // Finish current psi
            curPsi = dot(rDPtr[celli], curPsi);

            // Distribute the neighbour side using current psi
            for (label curFace=fStart; curFace<fEnd; curFace++)
            {
                bPrimePtr[uPtr[curFace]] -= dot(lowerPtr[curFace], curPsi);
            }

            psiPtr[celli] = curPsi;
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
