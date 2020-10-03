/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "freestreamPressureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    supersonic_(false)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    supersonic_
    (
        dict.lookupOrDefault<Switch>("supersonic", false)
    )
{
    freestreamValue() = scalarField("freestreamValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(freestreamValue());
    }

    refGrad() = Zero;
    valueFraction() = 0;
}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    UName_(psf.UName_),
    supersonic_(psf.supersonic_)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    UName_(psf.UName_),
    supersonic_(psf.supersonic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const Field<vector>& Up =
        patch().template lookupPatchField<volVectorField, vector>
        (
            UName_
        );

    const Field<scalar> magUp(mag(Up));

    const Field<vector> nf(patch().nf());

    Field<scalar>& vf = valueFraction();

    if (supersonic_)
    {
        forAll(vf, i)
        {
            if (magUp[i] > vSmall)
            {
                vf[i] = 0.5 - 0.5*(Up[i] & nf[i])/magUp[i];
            }
            else
            {
                vf[i] = 0.5;
            }
        }
    }
    else
    {
        forAll(vf, i)
        {
            if (magUp[i] > vSmall)
            {
                vf[i] = 0.5 + 0.5*(Up[i] & nf[i])/magUp[i];
            }
            else
            {
                vf[i] = 0.5;
            }
        }
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::freestreamPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry(os, "freestreamValue", freestreamValue());
    writeEntry(os, "supersonic", supersonic_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freestreamPressureFvPatchScalarField
    );
}

// ************************************************************************* //
