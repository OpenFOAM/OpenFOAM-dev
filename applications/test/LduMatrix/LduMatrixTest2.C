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

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "LduMatrix.H"
#include "tensorField.H"
#include "TPCG.H"
#include "TPBiCG.H"
#include "NoPreconditioner.H"

using namespace Foam;

typedef Foam::LduMatrix<vector, tensor, scalar>
    lduVectorMatrix;
defineNamedTemplateTypeNameAndDebug(lduVectorMatrix, 0);

typedef Foam::DiagonalSolver<vector, tensor, scalar>
    lduVectorDiagonalSolver;
defineNamedTemplateTypeNameAndDebug(lduVectorDiagonalSolver, 0);

template<>
const vector lduVectorMatrix::great_(1e15, 1e15, 1e15);

template<>
const vector lduVectorMatrix::small_(1e-15, 1e-15, 1e-15);

namespace Foam
{
    typedef LduMatrix<vector, tensor, scalar>::preconditioner
        lduVectorPreconditioner;
    defineTemplateRunTimeSelectionTable(lduVectorPreconditioner, symMatrix);
    defineTemplateRunTimeSelectionTable(lduVectorPreconditioner, asymMatrix);

    typedef LduMatrix<vector, tensor, scalar>::smoother
        lduVectorSmoother;
    defineTemplateRunTimeSelectionTable(lduVectorSmoother, symMatrix);
    defineTemplateRunTimeSelectionTable(lduVectorSmoother, asymMatrix);

    typedef LduMatrix<vector, tensor, scalar>::solver
        lduVectorSolver;
    defineTemplateRunTimeSelectionTable(lduVectorSolver, symMatrix);
    defineTemplateRunTimeSelectionTable(lduVectorSolver, asymMatrix);

    typedef TPCG<vector, tensor, scalar> TPCGVector;
    defineNamedTemplateTypeNameAndDebug(TPCGVector, 0);

    LduMatrix<vector, tensor, scalar>::solver::
        addsymMatrixConstructorToTable<TPCGVector>
        addTPCGSymMatrixConstructorToTable_;

    typedef TPBiCG<vector, tensor, scalar> TPBiCGVector;
    defineNamedTemplateTypeNameAndDebug(TPBiCGVector, 0);

    LduMatrix<vector, tensor, scalar>::solver::
        addasymMatrixConstructorToTable<TPBiCGVector>
        addTPBiCGSymMatrixConstructorToTable_;

    typedef NoPreconditioner<vector, tensor, scalar> NoPreconditionerVector;
    defineNamedTemplateTypeNameAndDebug(NoPreconditionerVector, 0);

    LduMatrix<vector, tensor, scalar>::preconditioner::
        addsymMatrixConstructorToTable<NoPreconditionerVector>
        addNoPreconditionerSymMatrixConstructorToTable_;

    LduMatrix<vector, tensor, scalar>::preconditioner::
        addasymMatrixConstructorToTable<NoPreconditionerVector>
        addNoPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    volVectorField psi
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    lduVectorMatrix testMatrix(mesh);
    testMatrix.diag() = 2*I;
    testMatrix.source() = pTraits<vector>::one;
    testMatrix.upper() = 0.1;
    testMatrix.lower() = -0.1;

    Info<< testMatrix << endl;

    FieldField<Field, scalar> boundaryCoeffs(0);
    FieldField<Field, scalar> internalCoeffs(0);

    autoPtr<lduVectorMatrix::solver> testMatrixSolver =
    lduVectorMatrix::solver::New
    (
        psi.name(),
        testMatrix,
        //boundaryCoeffs,
        //internalCoeffs,
        //psi.boundaryField().interfaces(),
        IStringStream
        (
            "PBiCG"
            "{"
            "    preconditioner   none;"
            "    tolerance        (1e-05 1e-05 1e-05);"
            "    relTol           (0 0 0);"
            "}"
        )()
    );

    lduVectorMatrix::solverPerformance solverPerf =
        testMatrixSolver->solve(psi);

    solverPerf.print();

    Info<< psi << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
