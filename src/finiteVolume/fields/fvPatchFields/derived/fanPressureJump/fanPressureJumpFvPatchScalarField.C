/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "fanPressureJumpFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fanPressureJumpFvPatchScalarField::calcFanJump()
{
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalar sign = reverse_ ? -1 : 1;

    if (fanCurve_.valid())
    {
        // Preferred method

        scalar volFlowRate = 0;

        if (phip.internalField().dimensions() == dimFlux)
        {
            volFlowRate = gSum(phip);
        }
        else if
        (
            phip.internalField().dimensions() == dimMassFlux
        )
        {
            const scalarField& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            volFlowRate = gSum(phip/rhop);
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of phi are not correct"
                << "\n    on patch " << patch().name()
                << " of field " << internalField().name()
                << " in file " << internalField().objectPath() << nl
                << exit(FatalError);
        }

        jump_ = sign*max(fanCurve_->value(max(sign*volFlowRate, 0)), 0);
    }
    else
    {
        // Backwards compatibility fallback

        scalarField Un(max(sign*phip/patch().magSf(), scalar(0)));

        if (phip.internalField().dimensions() == dimFlux)
        {
            // Do nothing
        }
        else if
        (
            phip.internalField().dimensions() == dimMassFlux
        )
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            Un /= rhop;
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of phi are not correct"
                << "\n    on patch " << patch().name()
                << " of field " << internalField().name()
                << " in file " << internalField().objectPath() << nl
                << exit(FatalError);
        }

        jump_ = sign*max(jumpTable_->value(Un), scalar(0));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureJumpFvPatchScalarField::fanPressureJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchScalarField(p, iF),
    fanCurve_(),
    jumpTable_(),
    reverse_(dict.lookupOrDefault<Switch>("reverse", false)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (cyclicPatch().owner())
    {
        if (dict.found("fanCurve"))
        {
            // Preferred entry name
            fanCurve_ = Function1<scalar>::New("fanCurve", dict);
        }
        else if (dict.found("jumpTable"))
        {
            // Backwards compatibility fallback
            jumpTable_ = Function1<scalar>::New("jumpTable", dict);
        }
        else
        {
            // Call function construction for error message
            fanCurve_ = Function1<scalar>::New("fanCurve", dict);
        }
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::fanPressureJumpFvPatchScalarField::fanPressureJumpFvPatchScalarField
(
    const fanPressureJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedJumpFvPatchScalarField(ptf, p, iF, mapper),
    fanCurve_(ptf.fanCurve_, false),
    jumpTable_(ptf.jumpTable_, false),
    reverse_(ptf.reverse_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::fanPressureJumpFvPatchScalarField::fanPressureJumpFvPatchScalarField
(
    const fanPressureJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchScalarField(ptf, iF),
    fanCurve_(ptf.fanCurve_, false),
    jumpTable_(ptf.jumpTable_, false),
    reverse_(ptf.reverse_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureJumpFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (cyclicPatch().owner())
    {
        calcFanJump();
    }

    fixedJumpFvPatchScalarField::updateCoeffs();
}


void Foam::fanPressureJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchScalarField::write(os);

    if (cyclicPatch().owner())
    {
        if (fanCurve_.valid())
        {
            writeEntry(os, fanCurve_());
        }
        else
        {
            writeEntry(os, jumpTable_());
        }
    }

    writeEntryIfDifferent<Switch>(os, "reverse", false, reverse_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureJumpFvPatchScalarField
    );
};


// ************************************************************************* //
