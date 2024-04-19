/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureJumpFvPatchScalarField::fanPressureJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchScalarField(p, iF, dict, true),
    phiName_
    (
        cyclicPatch().owner()
      ? dict.lookupOrDefault<word>("phi", "phi")
      : word::null
    ),
    rhoName_
    (
        cyclicPatch().owner()
      ? dict.lookupOrDefault<word>("rho", "rho")
      : word::null
    ),
    fanCurve_(),
    jumpTable_(),
    reverse_(dict.lookupOrDefault<Switch>("reverse", false))
{
    if (cyclicPatch().owner())
    {
        if (!dict.found("fanCurve") && dict.found("jumpTable"))
        {
            // Backwards compatibility fallback
            jumpTable_ =
                Function1<scalar>::New
                (
                    "jumpTable",
                    dimVelocity,
                    iF.dimensions(),
                    dict
                );
        }
        else
        {
            // Preferred entry name
            fanCurve_ =
                Function1<scalar>::New
                (
                    "fanCurve",
                    dimVolumetricFlux,
                    iF.dimensions(),
                    dict
                );
        }
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
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    fanCurve_(ptf.fanCurve_, false),
    jumpTable_(ptf.jumpTable_, false),
    reverse_(ptf.reverse_)
{}


Foam::fanPressureJumpFvPatchScalarField::fanPressureJumpFvPatchScalarField
(
    const fanPressureJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    fanCurve_(ptf.fanCurve_, false),
    jumpTable_(ptf.jumpTable_, false),
    reverse_(ptf.reverse_)
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
        const fvsPatchField<scalar>& phip =
            patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

        const scalar sign = reverse_ ? -1 : 1;

        if (fanCurve_.valid())
        {
            // Preferred method

            scalar volFlowRate = 0;

            if (phip.internalField().dimensions() == dimVolumetricFlux)
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

            jumpRef() = sign*max(fanCurve_->value(max(sign*volFlowRate, 0)), 0);
        }
        else
        {
            // Backwards compatibility fallback

            scalarField Un(max(sign*phip/patch().magSf(), scalar(0)));

            if (phip.internalField().dimensions() == dimVolumetricFlux)
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

            jumpRef() = sign*max(jumpTable_->value(Un), scalar(0));
        }
    }

    fixedJumpFvPatchScalarField::updateCoeffs();
}


void Foam::fanPressureJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchScalarField::write(os);

    if (cyclicPatch().owner())
    {
        writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

        if (fanCurve_.valid())
        {
            writeEntry(os, fanCurve_());
        }
        else
        {
            writeEntry(os, jumpTable_());
        }

        writeEntryIfDifferent<Switch>(os, "reverse", false, reverse_);
    }
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
