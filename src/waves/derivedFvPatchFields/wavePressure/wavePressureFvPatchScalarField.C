/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "wavePressureFvPatchScalarField.H"
#include "waveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "levelSet.H"
#include "volFields.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wavePressureFvPatchScalarField::wavePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::wavePressureFvPatchScalarField::wavePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::wavePressureFvPatchScalarField::wavePressureFvPatchScalarField
(
    const wavePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::wavePressureFvPatchScalarField::wavePressureFvPatchScalarField
(
    const wavePressureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::wavePressureFvPatchScalarField::wavePressureFvPatchScalarField
(
    const wavePressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::wavePressureFvPatchScalarField::p() const
{
    const scalar t = db().time().timeOutputValue();
    const waveSuperposition& waves = waveSuperposition::New(db());

    return
        levelSetAverage
        (
            patch(),
            waves.height(t, patch().Cf()),
            waves.height(t, patch().patch().localPoints()),
            waves.pGas(t, patch().Cf())(),
            waves.pGas(t, patch().patch().localPoints())(),
            waves.pLiquid(t, patch().Cf())(),
            waves.pLiquid(t, patch().patch().localPoints())()
        );
}


Foam::tmp<Foam::scalarField> Foam::wavePressureFvPatchScalarField::pn() const
{
    const scalar t = db().time().timeOutputValue();
    const waveSuperposition& waves = waveSuperposition::New(db());
    const waveVelocityFvPatchVectorField& Up =
        refCast<const waveVelocityFvPatchVectorField>
        (
            patch().lookupPatchField<volVectorField, scalar>(UName_)
        );

    const fvMeshSubset& subset = Up.faceCellSubset();
    const fvMesh& meshs = subset.subMesh();
    const label patchis = findIndex(subset.patchMap(), patch().index());

    const scalarField ps
    (
        levelSetAverage
        (
            meshs,
            waves.height(t, meshs.cellCentres())(),
            waves.height(t, meshs.points())(),
            waves.pGas(t, meshs.cellCentres())(),
            waves.pGas(t, meshs.points())(),
            waves.pLiquid(t, meshs.cellCentres())(),
            waves.pLiquid(t, meshs.points())()
        )
    );

    tmp<scalarField> tResult(new scalarField(patch().size()));
    scalarField& result = tResult.ref();

    if (patchis != -1)
    {
        forAll(meshs.boundary()[patchis], is)
        {
            const label fs = is + meshs.boundary()[patchis].patch().start();
            const label cs = meshs.boundary()[patchis].faceCells()[is];
            const label f = subset.faceMap()[fs];
            const label i = patch().patch().whichFace(f);
            result[i] = ps[cs];
        }
    }

    return tResult;
}


void Foam::wavePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, scalar>(UName_);

    if (!isA<waveVelocityFvPatchVectorField>(Up))
    {
        FatalErrorInFunction
            << "The corresponding condition for the velocity "
            << "field " << UName_ << " on patch " << patch().name()
            << " is not of type " << waveVelocityFvPatchVectorField::typeName
            << exit(FatalError);
    }

    const waveVelocityFvPatchVectorField& Uwp =
        refCast<const waveVelocityFvPatchVectorField>(Up);

    if (Uwp.pName() != internalField().name())
    {
        FatalErrorInFunction
            << "The corresponding condition for the velocity "
            << "field " << UName_ << " on patch " << patch().name()
            << " does not have the pressure set to " << internalField().name()
            << exit(FatalError);
    }

    const scalarField p(this->p()), pn(this->pn());
    const scalarField out(pos0(Uwp.U() & patch().Sf()));

    valueFraction() = out;
    refValue() = p;
    refGrad() = (p - pn)*patch().deltaCoeffs();

    if (internalField().dimensions() == dimPressure)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        refValue() *= rhop;
        refGrad() *= rhop;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::wavePressureFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, wavePressureFvPatchScalarField);
}

// ************************************************************************* //
