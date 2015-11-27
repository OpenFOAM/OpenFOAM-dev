/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF),
    relax_(1.0),
    fixedDmdt_(0.0),
    L_(0.0)
{
    checkType();
}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    fixedDmdt_(dict.lookupOrDefault<scalar>("fixedDmdt", 0.0)),
    L_(dict.lookupOrDefault<scalar>("L", 0.0))
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    fixedDmdt_(psf.fixedDmdt_),
    L_(psf.L_)
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    relax_(psf.relax_),
    fixedDmdt_(psf.fixedDmdt_),
    L_(psf.L_)
{}

alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    relax_(psf.relax_),
    fixedDmdt_(psf.fixedDmdt_),
    L_(psf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    dmdt_ = (1 - relax_)*dmdt_ + relax_*fixedDmdt_;
    mDotL_ = dmdt_*L_;

    operator==(calcAlphat(*this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;
    os.writeKeyword("fixedDmdt") << fixedDmdt_ << token::END_STATEMENT << nl;
    os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
    dmdt_.writeEntry("dmdt", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
