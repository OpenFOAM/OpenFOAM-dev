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

#include "LUscalarMatrix.H"
#include "lduMatrix.H"
#include "procLduMatrix.H"
#include "procLduInterface.H"
#include "cyclicLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LUscalarMatrix, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LUscalarMatrix::LUscalarMatrix()
:
    comm_(Pstream::worldComm)
{}


Foam::LUscalarMatrix::LUscalarMatrix(const scalarSquareMatrix& matrix)
:
    scalarSquareMatrix(matrix),
    comm_(Pstream::worldComm),
    pivotIndices_(m())
{
    LUDecompose(*this, pivotIndices_);
}


Foam::LUscalarMatrix::LUscalarMatrix
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    comm_(ldum.mesh().comm())
{
    if (Pstream::parRun())
    {
        PtrList<procLduMatrix> lduMatrices(Pstream::nProcs(comm_));

        label lduMatrixi = 0;

        lduMatrices.set
        (
            lduMatrixi++,
            new procLduMatrix
            (
                ldum,
                interfaceCoeffs,
                interfaces
            )
        );

        if (Pstream::master(comm_))
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave(comm_);
                slave++
            )
            {
                lduMatrices.set
                (
                    lduMatrixi++,
                    new procLduMatrix
                    (
                        IPstream
                        (
                            Pstream::commsTypes::scheduled,
                            slave,
                            0,          // bufSize
                            Pstream::msgType(),
                            comm_
                        )()
                    )
                );
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                0,              // bufSize
                Pstream::msgType(),
                comm_
            );
            procLduMatrix cldum
            (
                ldum,
                interfaceCoeffs,
                interfaces
            );
            toMaster<< cldum;

        }

        if (Pstream::master(comm_))
        {
            label nCells = 0;
            forAll(lduMatrices, i)
            {
                nCells += lduMatrices[i].size();
            }

            scalarSquareMatrix m(nCells, 0.0);
            transfer(m);
            convert(lduMatrices);
        }
    }
    else
    {
        label nCells = ldum.lduAddr().size();
        scalarSquareMatrix m(nCells, 0.0);
        transfer(m);
        convert(ldum, interfaceCoeffs, interfaces);
    }

    if (Pstream::master(comm_))
    {
        label mRows = m();
        label nColumns = n();

        if (debug)
        {
            Pout<< "LUscalarMatrix : size:" << mRows << endl;
            for (label rowI = 0; rowI < mRows; rowI++)
            {
                const scalar* row = operator[](rowI);

                Pout<< "cell:" << rowI << " diagCoeff:" << row[rowI] << endl;

                Pout<< "    connects to upper cells :";
                for (label columnI = rowI+1; columnI < nColumns; columnI++)
                {
                    if (mag(row[columnI]) > small)
                    {
                        Pout<< ' ' << columnI << " (coeff:" << row[columnI]
                            << ")";
                    }
                }
                Pout<< endl;
                Pout<< "    connects to lower cells :";
                for (label columnI = 0; columnI < rowI; columnI++)
                {
                    if (mag(row[columnI]) > small)
                    {
                        Pout<< ' ' << columnI << " (coeff:" << row[columnI]
                            << ")";
                    }
                }
                Pout<< endl;
            }
            Pout<< endl;
        }

        pivotIndices_.setSize(m());
        LUDecompose(*this, pivotIndices_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LUscalarMatrix::convert
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
{
    const label* __restrict__ uPtr = ldum.lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = ldum.lduAddr().lowerAddr().begin();

    const scalar* __restrict__ diagPtr = ldum.diag().begin();
    const scalar* __restrict__ upperPtr = ldum.upper().begin();
    const scalar* __restrict__ lowerPtr = ldum.lower().begin();

    const label nCells = ldum.diag().size();
    const label nFaces = ldum.upper().size();

    for (label cell=0; cell<nCells; cell++)
    {
        operator[](cell)[cell] = diagPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        label uCell = uPtr[face];
        label lCell = lPtr[face];

        operator[](uCell)[lCell] = lowerPtr[face];
        operator[](lCell)[uCell] = upperPtr[face];
    }

    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const lduInterface& interface = interfaces[inti].interface();

            // Assume any interfaces are cyclic ones

            const label* __restrict__ lPtr = interface.faceCells().begin();

            const cyclicLduInterface& cycInterface =
                refCast<const cyclicLduInterface>(interface);
            label nbrInt = cycInterface.neighbPatchID();
            const label* __restrict__ uPtr =
                interfaces[nbrInt].interface().faceCells().begin();

            const scalar* __restrict__ nbrUpperLowerPtr =
                interfaceCoeffs[nbrInt].begin();

            label inFaces = interface.faceCells().size();

            for (label face=0; face<inFaces; face++)
            {
                label uCell = lPtr[face];
                label lCell = uPtr[face];

                operator[](uCell)[lCell] -= nbrUpperLowerPtr[face];
            }
        }
    }
}


void Foam::LUscalarMatrix::convert
(
    const PtrList<procLduMatrix>& lduMatrices
)
{
    procOffsets_.setSize(lduMatrices.size() + 1);
    procOffsets_[0] = 0;

    forAll(lduMatrices, ldumi)
    {
        procOffsets_[ldumi+1] = procOffsets_[ldumi] + lduMatrices[ldumi].size();
    }

    forAll(lduMatrices, ldumi)
    {
        const procLduMatrix& lduMatrixi = lduMatrices[ldumi];
        label offset = procOffsets_[ldumi];

        const label* __restrict__ uPtr = lduMatrixi.upperAddr_.begin();
        const label* __restrict__ lPtr = lduMatrixi.lowerAddr_.begin();

        const scalar* __restrict__ diagPtr = lduMatrixi.diag_.begin();
        const scalar* __restrict__ upperPtr = lduMatrixi.upper_.begin();
        const scalar* __restrict__ lowerPtr = lduMatrixi.lower_.begin();

        const label nCells = lduMatrixi.size();
        const label nFaces = lduMatrixi.upper_.size();

        for (label cell=0; cell<nCells; cell++)
        {
            label globalCell = cell + offset;
            operator[](globalCell)[globalCell] = diagPtr[cell];
        }

        for (label face=0; face<nFaces; face++)
        {
            label uCell = uPtr[face] + offset;
            label lCell = lPtr[face] + offset;

            operator[](uCell)[lCell] = lowerPtr[face];
            operator[](lCell)[uCell] = upperPtr[face];
        }

        const PtrList<procLduInterface>& interfaces =
            lduMatrixi.interfaces_;

        forAll(interfaces, inti)
        {
            const procLduInterface& interface = interfaces[inti];

            if (interface.myProcNo_ == interface.neighbProcNo_)
            {
                const label* __restrict__ ulPtr = interface.faceCells_.begin();

                const scalar* __restrict__ upperLowerPtr =
                    interface.coeffs_.begin();

                label inFaces = interface.faceCells_.size()/2;

                for (label face=0; face<inFaces; face++)
                {
                    label uCell = ulPtr[face] + offset;
                    label lCell = ulPtr[face + inFaces] + offset;

                    operator[](uCell)[lCell] -= upperLowerPtr[face + inFaces];
                    operator[](lCell)[uCell] -= upperLowerPtr[face];
                }
            }
            else if (interface.myProcNo_ < interface.neighbProcNo_)
            {
                // Interface to neighbour proc. Find on neighbour proc the
                // corresponding interface. The problem is that there can
                // be multiple interfaces between two processors (from
                // processorCyclics) so also compare the communication tag

                const PtrList<procLduInterface>& neiInterfaces =
                    lduMatrices[interface.neighbProcNo_].interfaces_;

                label neiInterfacei = -1;

                forAll(neiInterfaces, ninti)
                {
                    if
                    (
                        (
                            neiInterfaces[ninti].neighbProcNo_
                         == interface.myProcNo_
                        )
                     && (neiInterfaces[ninti].tag_ ==  interface.tag_)
                    )
                    {
                        neiInterfacei = ninti;
                        break;
                    }
                }

                if (neiInterfacei == -1)
                {
                    FatalErrorInFunction << exit(FatalError);
                }

                const procLduInterface& neiInterface =
                    neiInterfaces[neiInterfacei];

                const label* __restrict__ uPtr = interface.faceCells_.begin();
                const label* __restrict__ lPtr =
                    neiInterface.faceCells_.begin();

                const scalar* __restrict__ upperPtr = interface.coeffs_.begin();
                const scalar* __restrict__ lowerPtr =
                    neiInterface.coeffs_.begin();

                label inFaces = interface.faceCells_.size();
                label neiOffset = procOffsets_[interface.neighbProcNo_];

                for (label face=0; face<inFaces; face++)
                {
                    label uCell = uPtr[face] + offset;
                    label lCell = lPtr[face] + neiOffset;

                    operator[](uCell)[lCell] -= lowerPtr[face];
                    operator[](lCell)[uCell] -= upperPtr[face];
                }
            }
        }
    }
}


void Foam::LUscalarMatrix::printDiagonalDominance() const
{
    for (label i=0; i<m(); i++)
    {
        scalar sum = 0.0;
        for (label j=0; j<m(); j++)
        {
            if (i != j)
            {
                sum += operator[](i)[j];
            }
        }
        Info<< mag(sum)/mag(operator[](i)[i]) << endl;
    }
}


void Foam::LUscalarMatrix::decompose(const scalarSquareMatrix& M)
{
    scalarSquareMatrix::operator=(M);
    pivotIndices_.setSize(m());
    LUDecompose(*this, pivotIndices_);
}


void Foam::LUscalarMatrix::inv(scalarSquareMatrix& M) const
{
    scalarField source(m());

    for (label j=0; j<m(); j++)
    {
        source = Zero;
        source[j] = 1;
        LUBacksubstitute(*this, pivotIndices_, source);
        for (label i=0; i<m(); i++)
        {
            M(i, j) = source[i];
        }
    }
}


// ************************************************************************* //
