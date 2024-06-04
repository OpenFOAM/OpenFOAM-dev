/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "epsilonmWallFunctionFvPatchScalarField.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::epsilonmWallFunctionFvPatchScalarField::manipulateMatrixMaster
(
    fvMatrix<scalar>& matrix
)
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    matrix.setValues
    (
        wallCells(),
        UIndirectList<scalar>(internalField(), wallCells()),
        wallCellFraction()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonmWallFunctionFvPatchScalarField::
epsilonmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallCellWallFunctionFvPatchScalarField(p, iF, dict)
{}


Foam::epsilonmWallFunctionFvPatchScalarField::
epsilonmWallFunctionFvPatchScalarField
(
    const epsilonmWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    wallCellWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::epsilonmWallFunctionFvPatchScalarField::
epsilonmWallFunctionFvPatchScalarField
(
    const epsilonmWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallCellWallFunctionFvPatchScalarField(ewfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::epsilonmWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    initMaster();
}


void Foam::epsilonmWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    if (masterPatchIndex() == -1)
    {
        FatalErrorInFunction
            << "updateCoeffs must be called before manipulateMatrix"
            << exit(FatalError);
    }

    manipulateMatrixMaster(matrix);

    fvPatchScalarField::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonmWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
