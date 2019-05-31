/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    vaporPhaseName_("vapor"),
    relax_(1.0),
    fixedDmdt_(0.0),
    L_(0.0)
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    vaporPhaseName_(dict.lookup("vaporPhase")),
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

bool alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
activePhasePair(const phasePairKey& phasePair) const
{
    if (phasePair == phasePairKey(vaporPhaseName_, internalField().group()))
    {
        return true;
    }
    else
    {
        return false;
    }
}


const scalarField& alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
dmdt(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return dmdt_;
    }
    else
    {
        FatalErrorInFunction
            << " dmdt requested for invalid phasePair!"
            << abort(FatalError);

        return mDotL_;
    }
}


const scalarField& alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
mDotL(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return mDotL_;
    }
    else
    {
        FatalErrorInFunction
            << " mDotL requested for invalid phasePair!"
            << abort(FatalError);

        return mDotL_;
    }
}


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
    writeEntry(os, "vaporPhase", vaporPhaseName_);
    writeEntry(os, "relax", relax_);
    writeEntry(os, "fixedDmdt", fixedDmdt_);
    writeEntry(os, "L", L_);
    writeEntry(os, "dmdt", dmdt_);
    writeEntry(os, "value", *this);
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
