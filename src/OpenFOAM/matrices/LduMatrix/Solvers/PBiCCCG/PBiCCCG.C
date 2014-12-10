/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "PBiCCCG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::PBiCCCG<Type, DType, LUType>::PBiCCCG
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
Foam::SolverPerformance<Type>
Foam::PBiCCCG<Type, DType, LUType>::solve
(
    Field<Type>& psi
) const
{
    word preconditionerName(this->controlDict_.lookup("preconditioner"));

    // --- Setup class containing solver performance data
    SolverPerformance<Type> solverPerf
    (
        preconditionerName + typeName,
        this->fieldName_
    );

    register label nCells = psi.size();

    Type* __restrict__ psiPtr = psi.begin();

    Field<Type> pA(nCells);
    Type* __restrict__ pAPtr = pA.begin();

    Field<Type> pT(nCells, pTraits<Type>::zero);
    Type* __restrict__ pTPtr = pT.begin();

    Field<Type> wA(nCells);
    Type* __restrict__ wAPtr = wA.begin();

    Field<Type> wT(nCells);
    Type* __restrict__ wTPtr = wT.begin();

    scalar wArT = 1e15; //this->matrix_.great_;
    scalar wArTold = wArT;

    // --- Calculate A.psi and T.psi
    this->matrix_.Amul(wA, psi);
    this->matrix_.Tmul(wT, psi);

    // --- Calculate initial residual and transpose residual fields
    Field<Type> rA(this->matrix_.source() - wA);
    Field<Type> rT(this->matrix_.source() - wT);
    Type* __restrict__ rAPtr = rA.begin();
    Type* __restrict__ rTPtr = rT.begin();

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
            // --- Store previous wArT
            wArTold = wArT;

            // --- Precondition residuals
            preconPtr->precondition(wA, rA);
            preconPtr->preconditionT(wT, rT);

            // --- Update search directions:
            wArT = gSumProd(wA, rT);

            if (solverPerf.nIterations() == 0)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                    pTPtr[cell] = wTPtr[cell];
                }
            }
            else
            {
                scalar beta = wArT/wArTold;

                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + (beta* pAPtr[cell]);
                    pTPtr[cell] = wTPtr[cell] + (beta* pTPtr[cell]);
                }
            }


            // --- Update preconditioned residuals
            this->matrix_.Amul(wA, pA);
            this->matrix_.Tmul(wT, pT);

            scalar wApT = gSumProd(wA, pT);

            // --- Test for singularity
            if
            (
                solverPerf.checkSingularity
                (
                    cmptDivide(pTraits<Type>::one*mag(wApT), normFactor)
                )
            )
            {
                break;
            }


            // --- Update solution and residual:

            scalar alpha = wArT/wApT;

            for (register label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += (alpha* pAPtr[cell]);
                rAPtr[cell] -= (alpha* wAPtr[cell]);
                rTPtr[cell] -= (alpha* wTPtr[cell]);
            }

            solverPerf.finalResidual() =
                cmptDivide(gSumCmptMag(rA), normFactor);

        } while
        (
            (
                solverPerf.nIterations()++ < this->maxIter_
            && !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
            )
         || solverPerf.nIterations() < this->minIter_
        );
    }

    return solverPerf;
}


// ************************************************************************* //
