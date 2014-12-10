/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "porousBafflePressureFvPatchField.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    D_(0),
    I_(0),
    length_(0)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    D_(readScalar(dict.lookup("D"))),
    I_(readScalar(dict.lookup("I"))),
    length_(readScalar(dict.lookup("length")))
{
    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<scalar>(ptf),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousBafflePressureFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>("phi");

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    scalarField Un(phip/patch().magSf());

    scalarField magUn(mag(Un));

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    jump_ =
        -sign(Un)
        *(
            I_*turbModel.nuEff(patch().index())
          + D_*0.5*magUn
         )*magUn*length_;

    if (debug)
    {
        scalar avePressureJump = gAverage(jump_);
        scalar aveVelocity = gAverage(mag(Un));

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << " Average pressure drop :" << avePressureJump
            << " Average velocity :" << aveVelocity
            << endl;
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::porousBafflePressureFvPatchField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    os.writeKeyword("D") << D_ << token::END_STATEMENT << nl;
    os.writeKeyword("I") << I_ << token::END_STATEMENT << nl;
    os.writeKeyword("length") << length_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousBafflePressureFvPatchField
    );
}

// ************************************************************************* //
