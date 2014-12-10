/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "GAMGSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::interpolate
(
    scalarField& psi,
    scalarField& Apsi,
    const lduMatrix& m,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    scalar* __restrict__ psiPtr = psi.begin();

    const label* const __restrict__ uPtr = m.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = m.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ diagPtr = m.diag().begin();
    const scalar* const __restrict__ upperPtr = m.upper().begin();
    const scalar* const __restrict__ lowerPtr = m.lower().begin();

    Apsi = 0;
    scalar* __restrict__ ApsiPtr = Apsi.begin();

    m.initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    register const label nFaces = m.upper().size();
    for (register label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }

    m.updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    register const label nCells = m.diag().size();
    for (register label celli=0; celli<nCells; celli++)
    {
        psiPtr[celli] = -ApsiPtr[celli]/(diagPtr[celli]);
    }
}


void Foam::GAMGSolver::interpolate
(
    scalarField& psi,
    scalarField& Apsi,
    const lduMatrix& m,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const labelList& restrictAddressing,
    const scalarField& psiC,
    const direction cmpt
) const
{
    interpolate
    (
        psi,
        Apsi,
        m,
        interfaceBouCoeffs,
        interfaces,
        cmpt
    );

    register const label nCells = m.diag().size();
    scalar* __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = m.diag().begin();

    register const label nCCells = psiC.size();
    scalarField corrC(nCCells, 0);
    scalarField diagC(nCCells, 0);

    for (register label celli=0; celli<nCells; celli++)
    {
        corrC[restrictAddressing[celli]] += diagPtr[celli]*psiPtr[celli];
        diagC[restrictAddressing[celli]] += diagPtr[celli];
    }

    for (register label ccelli=0; ccelli<nCCells; ccelli++)
    {
        corrC[ccelli] = psiC[ccelli] - corrC[ccelli]/diagC[ccelli];
    }

    for (register label celli=0; celli<nCells; celli++)
    {
        psiPtr[celli] += corrC[restrictAddressing[celli]];
    }
}


// ************************************************************************* //
