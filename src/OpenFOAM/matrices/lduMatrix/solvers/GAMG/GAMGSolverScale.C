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

#include "GAMGSolver.H"
#include "vector2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::scale
(
    scalarField& field,
    scalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const scalarField& source,
    const direction cmpt
) const
{
    A.Amul
    (
        Acf,
        field,
        interfaceLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );

    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    forAll(field, i)
    {
        scalingFactorNum += source[i]*field[i];
        scalingFactorDenom += Acf[i]*field[i];
    }

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    A.mesh().reduce(scalingVector, sumOp<vector2D>());

    scalar sf = scalingVector.x()/stabilise(scalingVector.y(), VSMALL);

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalarField& D = A.diag();

    forAll(field, i)
    {
        field[i] = sf*field[i] + (source[i] - sf*Acf[i])/D[i];
    }
}


// ************************************************************************* //
