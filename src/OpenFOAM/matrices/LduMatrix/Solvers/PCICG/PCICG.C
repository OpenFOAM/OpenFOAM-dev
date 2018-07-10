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

#include "PCICG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::PCICG<Type, DType, LUType>::PCICG
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    LduMatrix<Type, DType, LUType>::solver
    (
        fieldName,
        matrix,
        solverDict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
typename Foam::SolverPerformance<Type>
Foam::PCICG<Type, DType, LUType>::solve(Field<Type>& psi) const
{
    word preconditionerName(this->controlDict_.lookup("preconditioner"));

    // --- Setup class containing solver performance data
    SolverPerformance<Type> solverPerf
    (
        preconditionerName + typeName,
        this->fieldName_
    );

    label nIter = 0;

    label nCells = psi.size();

    Type* __restrict__ psiPtr = psi.begin();

    Field<Type> pA(nCells);
    Type* __restrict__ pAPtr = pA.begin();

    Field<Type> wA(nCells);
    Type* __restrict__ wAPtr = wA.begin();

    Type wArA = solverPerf.great_*pTraits<Type>::one;
    Type wArAold = wArA;

    // --- Calculate A.psi
    this->matrix_.Amul(wA, psi);

    // --- Calculate initial residual field
    Field<Type> rA(this->matrix_.source() - wA);
    Type* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    Type normFactor = this->normFactor(psi, wA, pA);

    if (LduMatrix<Type, DType, LUType>::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = cmptDivide(gSumCmptMag(rA), normFactor);
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        this->minIter_ > 0
     || !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
    )
    {
        // --- Select and construct the preconditioner
        autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        preconPtr = LduMatrix<Type, DType, LUType>::preconditioner::New
        (
            *this,
            this->controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            preconPtr->precondition(wA, rA);

            // --- Update search directions:
            wArA = gSumCmptProd(wA, rA);

            if (nIter == 0)
            {
                for (label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                }
            }
            else
            {
                Type beta = cmptDivide
                (
                    wArA,
                    stabilise(wArAold, solverPerf.vsmall_)
                );

                for (label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + cmptMultiply(beta, pAPtr[cell]);
                }
            }


            // --- Update preconditioned residual
            this->matrix_.Amul(wA, pA);

            Type wApA = gSumCmptProd(wA, pA);


            // --- Test for singularity
            if
            (
                solverPerf.checkSingularity
                (
                    cmptDivide(cmptMag(wApA), normFactor)
                )
            )
            {
                break;
            }


            // --- Update solution and residual:

            Type alpha = cmptDivide
            (
                wArA,
                stabilise(wApA, solverPerf.vsmall_)
            );

            for (label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += cmptMultiply(alpha, pAPtr[cell]);
                rAPtr[cell] -= cmptMultiply(alpha, wAPtr[cell]);
            }

            solverPerf.finalResidual() =
                cmptDivide(gSumCmptMag(rA), normFactor);

        } while
        (
            (
                nIter++ < this->maxIter_
            && !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
            )
         || nIter < this->minIter_
        );
    }

    solverPerf.nIterations() =
        pTraits<typename pTraits<Type>::labelType>::one*nIter;

    return solverPerf;
}


// ************************************************************************* //
