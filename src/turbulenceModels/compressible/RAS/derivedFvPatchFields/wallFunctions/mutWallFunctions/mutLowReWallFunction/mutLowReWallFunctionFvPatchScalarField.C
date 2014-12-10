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

#include "mutLowReWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutLowReWallFunctionFvPatchScalarField::calcMut() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const mutLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const mutLowReWallFunctionFvPatchScalarField& mlrwfpsf
)
:
    mutWallFunctionFvPatchScalarField(mlrwfpsf)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const mutLowReWallFunctionFvPatchScalarField& mlrwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(mlrwfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> mutLowReWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbModel.y()[patchi];
    const scalarField nuw
    (
        turbModel.mu().boundaryField()[patchi]
       /turbModel.rho().boundaryField()[patchi]
    );
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    return y*sqrt(nuw*mag(Uw.snGrad()))/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
