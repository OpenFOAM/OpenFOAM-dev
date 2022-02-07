/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "totalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pressureInletOutletVelocityFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicPressureFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi")
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    dynamicPressureFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    dynamicPressureFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_)
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicPressureFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalPressureFvPatchScalarField::updateCoeffs()
{
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    if (isA<pressureInletOutletVelocityFvPatchVectorField>(Up))
    {
        const pressureInletOutletVelocityFvPatchVectorField& Upiov =
            refCast<const pressureInletOutletVelocityFvPatchVectorField>(Up);

        if (Upiov.tangentialVelocity().valid())
        {
            const scalar t = this->db().time().userTimeValue();

            dynamicPressureFvPatchScalarField::updateCoeffs
            (
                p0_,
                0.5*neg(phip)*magSqr(Upiov.tangentialVelocity()->value(t))
              - 0.5*neg(phip)*magSqr(Up)
            );

            return;
        }
    }

    dynamicPressureFvPatchScalarField::updateCoeffs
    (
        p0_,
        -0.5*neg(phip)*magSqr(Up)
    );
}


void Foam::totalPressureFvPatchScalarField::write(Ostream& os) const
{
    dynamicPressureFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
