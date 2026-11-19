/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "pressureInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    tangentialVelocity_
    (
        iF.name(),
        "tangentialVelocity",
        p,
        iF.dimensions(),
        refValue(),
        dict,
        Zero
    ),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    const vectorField n(patch().nf());

    refValue() -= n*(n & refValue());
    refGrad() = Zero;
    valueFraction() = Zero;

    if (dict.found("value"))
    {
        vectorField::operator=
        (
            vectorField("value", iF.dimensions(), dict, p.size())
        );
    }
    else if (p.time().completeCase())
    {
        const scalarField phip(patchInternalField() & n);
        valueFraction() = neg(phip)*(I - sqr(n));

        directionMixedFvPatchVectorField::updateCoeffs();
        directionMixedFvPatchVectorField::evaluate();
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unable to evaluate function for incomplete case "
                "and 'value' entry not provided."
            << exit(FatalIOError);
    }
}


Foam::pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    tangentialVelocity_
    (
        ptf.tangentialVelocity_,
        p,
        refValue()
    ),
    phiName_(ptf.phiName_)
{}


Foam::pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    tangentialVelocity_
    (
        ptf.tangentialVelocity_,
        refValue()
    ),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureInletOutletVelocityFvPatchVectorField::map
(
    const fvPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    directionMixedFvPatchVectorField::map(ptf, mapper);
    tangentialVelocity_.map(!mapper.direct());
}


void Foam::pressureInletOutletVelocityFvPatchVectorField::reset
(
    const fvPatchField<vector>& ptf
)
{
    directionMixedFvPatchVectorField::reset(ptf);
    tangentialVelocity_.reset();
}


void Foam::pressureInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField n(patch().nf());

    if (tangentialVelocity_.update())
    {
        refValue() -= n*(n & refValue());
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    valueFraction() = neg(phip)*(I - sqr(n));

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::pressureInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
)
const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, tangentialVelocity_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureInletOutletVelocityFvPatchVectorField::operator=
(
    const expression::Assignment<vectorField>& a
)
{
    tmp<vectorField> normalValue = transform(valueFraction(), refValue());

    tmp<vectorField> pvf(new vectorField(this->size()));
    a.assign(pvf.ref());

    tmp<vectorField> transformGradValue = transform(I - valueFraction(), pvf);

    vectorField::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
