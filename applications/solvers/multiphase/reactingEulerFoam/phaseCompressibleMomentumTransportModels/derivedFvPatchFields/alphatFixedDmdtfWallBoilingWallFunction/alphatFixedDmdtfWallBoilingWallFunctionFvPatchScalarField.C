/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(p, iF),
    fixedDmdtf_(0),
    L_(0)
{}


alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(p, iF, dict),
    fixedDmdtf_(dict.lookupOrDefault<scalar>("fixedDmdtf", 0)),
    L_(dict.lookupOrDefault<scalar>("L", 0))
{}


alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    fixedDmdtf_(psf.fixedDmdtf_),
    L_(psf.L_)
{}


alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(psf),
    fixedDmdtf_(psf.fixedDmdtf_),
    L_(psf.L_)
{}


alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(psf, iF),
    fixedDmdtf_(psf.fixedDmdtf_),
    L_(psf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    dmdtf_ = (1 - relax_)*dmdtf_ + relax_*fixedDmdtf_;
    dmdtLf_ = dmdtf_*L_;

    operator==(calcAlphat(*this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    alphatPhaseChangeWallFunctionFvPatchScalarField::write(os);

    writeEntry(os, "fixedDmdtf", fixedDmdtf_);
    writeEntry(os, "L", L_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
