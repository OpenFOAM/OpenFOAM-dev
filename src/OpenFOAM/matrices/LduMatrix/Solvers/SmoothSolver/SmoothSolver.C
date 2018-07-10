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

#include "SmoothSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::SmoothSolver<Type, DType, LUType>::SmoothSolver
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
    ),
    nSweeps_(1)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::SmoothSolver<Type, DType, LUType>::readControls()
{
    LduMatrix<Type, DType, LUType>::solver::readControls();
    this->readControl(this->controlDict_, nSweeps_, "nSweeps");
}


template<class Type, class DType, class LUType>
Foam::SolverPerformance<Type>
Foam::SmoothSolver<Type, DType, LUType>::solve(Field<Type>& psi) const
{
    // --- Setup class containing solver performance data
    SolverPerformance<Type> solverPerf
    (
        typeName,
        this->fieldName_
    );

    label nIter = 0;

    // If the nSweeps_ is negative do a fixed number of sweeps
    if (nSweeps_ < 0)
    {
        autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        smootherPtr = LduMatrix<Type, DType, LUType>::smoother::New
        (
            this->fieldName_,
            this->matrix_,
            this->controlDict_
        );

        smootherPtr->smooth(psi, -nSweeps_);

        nIter -= nSweeps_;
    }
    else
    {
        Type normFactor = Zero;

        {
            Field<Type> Apsi(psi.size());
            Field<Type> temp(psi.size());

            // Calculate A.psi
            this->matrix_.Amul(Apsi, psi);

            // Calculate normalisation factor
            normFactor = this->normFactor(psi, Apsi, temp);

            // Calculate residual magnitude
            solverPerf.initialResidual() = cmptDivide
            (
                gSumCmptMag(this->matrix_.source() - Apsi),
                normFactor
            );
            solverPerf.finalResidual() = solverPerf.initialResidual();
        }

        if (LduMatrix<Type, DType, LUType>::debug >= 2)
        {
            Info<< "   Normalisation factor = " << normFactor << endl;
        }


        // Check convergence, solve if not converged
        if
        (
            this->minIter_ > 0
         || !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
        )
        {
            autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
            smootherPtr = LduMatrix<Type, DType, LUType>::smoother::New
            (
                this->fieldName_,
                this->matrix_,
                this->controlDict_
            );

            // Smoothing loop
            do
            {
                smootherPtr->smooth
                (
                    psi,
                    nSweeps_
                );

                // Calculate the residual to check convergence
                solverPerf.finalResidual() = cmptDivide
                (
                    gSumCmptMag(this->matrix_.residual(psi)),
                    normFactor
                );
            } while
            (
                (
                    (nIter += nSweeps_) < this->maxIter_
                && !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
                )
             || nIter < this->minIter_
            );
        }
    }

    solverPerf.nIterations() =
        pTraits<typename pTraits<Type>::labelType>::one*nIter;

    return solverPerf;
}


// ************************************************************************* //
