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

#include "LduMatrix.H"
#include "DiagonalSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::solver>
Foam::LduMatrix<Type, DType, LUType>::solver::New
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
{
    word solverName = solverDict.lookup("solver");

    if (matrix.diagonal())
    {
        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            new DiagonalSolver<Type, DType, LUType>
            (
                fieldName,
                matrix,
                solverDict
            )
        );
    }
    else if (matrix.symmetric())
    {
        typename symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(solverName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "LduMatrix<Type, DType, LUType>::solver::New", solverDict
            )   << "Unknown symmetric matrix solver " << solverName
                << endl << endl
                << "Valid symmetric matrix solvers are :" << endl
                << symMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                solverDict
            )
        );
    }
    else if (matrix.asymmetric())
    {
        typename asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(solverName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "LduMatrix<Type, DType, LUType>::solver::New", solverDict
            )   << "Unknown asymmetric matrix solver " << solverName
                << endl << endl
                << "Valid asymmetric matrix solvers are :" << endl
                << asymMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                solverDict
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "LduMatrix<Type, DType, LUType>::solver::New", solverDict
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::solver::solver
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    fieldName_(fieldName),
    matrix_(matrix),

    controlDict_(solverDict),

    maxIter_(1000),
    minIter_(0),
    tolerance_(1e-6*pTraits<Type>::one),
    relTol_(pTraits<Type>::zero)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::solver::readControls()
{
    readControl(controlDict_, maxIter_, "maxIter");
    readControl(controlDict_, minIter_, "minIter");
    readControl(controlDict_, tolerance_, "tolerance");
    readControl(controlDict_, relTol_, "relTol");
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::solver::read
(
    const dictionary& solverDict
)
{
    controlDict_ = solverDict;
    readControls();
}


template<class Type, class DType, class LUType>
Type Foam::LduMatrix<Type, DType, LUType>::solver::normFactor
(
    const Field<Type>& psi,
    const Field<Type>& Apsi,
    Field<Type>& tmpField
) const
{
    // --- Calculate A dot reference value of psi
    matrix_.sumA(tmpField);
    cmptMultiply(tmpField, tmpField, gAverage(psi));

    return stabilise
    (
        gSum(cmptMag(Apsi - tmpField) + cmptMag(matrix_.source() - tmpField)),
        SolverPerformance<Type>::small_
    );

    // At convergence this simpler method is equivalent to the above
    // return stabilise(2*gSumCmptMag(matrix_.source()), matrix_.small_);
}


// ************************************************************************* //
