/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "fixedNormalInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedNormalInletOutletVelocityFvPatchVectorField::
fixedNormalInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    fixTangentialInflow_(dict.lookup("fixTangentialInflow")),
    normalVelocity_
    (
        fvPatchVectorField::New(p, iF, dict.subDict("normalVelocity"))
    )
{
    fvPatchVectorField::operator=
    (
        vectorField("value", iF.dimensions(), dict, p.size())
    );
    refValue() = normalVelocity();
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::fixedNormalInletOutletVelocityFvPatchVectorField::
fixedNormalInletOutletVelocityFvPatchVectorField
(
    const fixedNormalInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    fixTangentialInflow_(ptf.fixTangentialInflow_),
    normalVelocity_
    (
        fvPatchVectorField::New(ptf.normalVelocity(), p, iF, mapper)
    )
{}


Foam::fixedNormalInletOutletVelocityFvPatchVectorField::
fixedNormalInletOutletVelocityFvPatchVectorField
(
    const fixedNormalInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    fixTangentialInflow_(pivpvf.fixTangentialInflow_),
    normalVelocity_(pivpvf.normalVelocity().clone(iF))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedNormalInletOutletVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    directionMixedFvPatchVectorField::map(ptf, mapper);

    const fixedNormalInletOutletVelocityFvPatchVectorField& fniovptf =
        refCast<const fixedNormalInletOutletVelocityFvPatchVectorField>(ptf);

    mapper(normalVelocity_.ref(), fniovptf.normalVelocity());
}


void Foam::fixedNormalInletOutletVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    directionMixedFvPatchVectorField::reset(ptf);

    const fixedNormalInletOutletVelocityFvPatchVectorField& fniovptf =
        refCast<const fixedNormalInletOutletVelocityFvPatchVectorField>(ptf);

    normalVelocity_->reset(fniovptf.normalVelocity());
}


void Foam::fixedNormalInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    normalVelocity_->evaluate();
    refValue() = normalVelocity();

    valueFraction() = sqr(patch().nf());

    if (fixTangentialInflow_)
    {
        const fvsPatchField<scalar>& phip =
            patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

        valueFraction() += neg(phip)*(I - valueFraction());
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::fixedNormalInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
)
const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "fixTangentialInflow", fixTangentialInflow_);
    writeKeyword(os, "normalVelocity")
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;
    normalVelocity_->write(os);
    os << decrIndent << indent << token::END_BLOCK << endl;
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fixedNormalInletOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    tmp<vectorField> normalValue = transform(valueFraction(), refValue());
    tmp<vectorField> transformGradValue = transform(I - valueFraction(), pvf);
    fvPatchField<vector>::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fixedNormalInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
