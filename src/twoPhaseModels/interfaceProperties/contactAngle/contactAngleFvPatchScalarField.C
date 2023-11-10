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

#include "contactAngleFvPatchScalarField.H"
#include "fieldMapper.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::contactAngleFvPatchScalarField::limitControls,
        4
    >::names[] =
    {
        "none",
        "gradient",
        "zeroGradient",
        "alpha"
    };
}


const Foam::NamedEnum
<
    Foam::contactAngleFvPatchScalarField::limitControls,
    4
> Foam::contactAngleFvPatchScalarField::limitControlNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contactAngleFvPatchScalarField::contactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    contactAngle_(contactAngleModel::New(dict)),
    limit_(limitControlNames_.read(dict.lookup("limit")))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::contactAngleFvPatchScalarField::contactAngleFvPatchScalarField
(
    const contactAngleFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(acpsf, p, iF, mapper),
    contactAngle_(acpsf.contactAngle_, false),
    limit_(acpsf.limit_)
{}


Foam::contactAngleFvPatchScalarField::contactAngleFvPatchScalarField
(
    const contactAngleFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(acpsf, iF),
    contactAngle_(acpsf.contactAngle_, false),
    limit_(acpsf.limit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::contactAngleFvPatchScalarField::cosTheta
(
    const fvPatchVectorField& Up,
    const vectorField& nHat
) const
{
    return contactAngle_->cosTheta(Up, nHat);
}


void Foam::contactAngleFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (limit_ == lcGradient)
    {
        gradient() =
        patch().deltaCoeffs()
       *(
           max(min
           (
               *this + gradient()/patch().deltaCoeffs(),
               scalar(1)), scalar(0)
           ) - *this
       );
    }
    else if (limit_ == lcZeroGradient)
    {
        gradient() = 0.0;
    }

    fixedGradientFvPatchScalarField::evaluate();

    if (limit_ == lcAlpha)
    {
        scalarField::operator=(max(min(*this, scalar(1)), scalar(0)));
    }
}


void Foam::contactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, contactAngle_());
    writeEntry(os, "limit", limitControlNames_[limit_]);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        contactAngleFvPatchScalarField
    );
}


// ************************************************************************* //
