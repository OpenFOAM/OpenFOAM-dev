/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "semiPermeableBaffleMassFractionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        semiPermeableBaffleMassFractionFvPatchScalarField::input,
        4
    >::names[] =
    {
        "none",
        "massFraction",
        "moleFraction",
        "partialPressure",
    };
}

const Foam::NamedEnum
<
    Foam::semiPermeableBaffleMassFractionFvPatchScalarField::input,
    4
> Foam::semiPermeableBaffleMassFractionFvPatchScalarField::inputNames_;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::basicSpecieMixture&
Foam::semiPermeableBaffleMassFractionFvPatchScalarField::composition
(
    const objectRegistry& db
)
{
    const word& name = basicThermo::dictName;

    if (db.foundObject<psiReactionThermo>(name))
    {
        return db.lookupObject<psiReactionThermo>(name).composition();
    }
    else if (db.foundObject<rhoReactionThermo>(name))
    {
        return db.lookupObject<rhoReactionThermo>(name).composition();
    }
    else
    {
        FatalErrorInFunction
            << "Could not find a multi-component thermodynamic model."
            << exit(FatalError);

        return NullObjectRef<basicSpecieMixture>();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    c_(0),
    input_(none),
    phiName_("phi"),
    pName_("p")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), NEARESTPATCHFACE, dict),
    mixedFvPatchScalarField(p, iF),
    c_(dict.lookupOrDefault<scalar>("c", scalar(0))),
    input_
    (
        c_ == scalar(0)
      ? none
      : inputNames_.read(dict.lookup("input"))
    ),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    pName_(dict.lookupOrDefault<word>("p", "p"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase(p.patch(), ptf),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    c_(ptf.c_),
    input_(ptf.input_),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_)
{}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf),
    c_(ptf.c_),
    input_(ptf.input_),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_)
{}


Foam::semiPermeableBaffleMassFractionFvPatchScalarField::
semiPermeableBaffleMassFractionFvPatchScalarField
(
    const semiPermeableBaffleMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    c_(ptf.c_),
    input_(ptf.input_),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::semiPermeableBaffleMassFractionFvPatchScalarField::phiY() const
{
    if (c_ == scalar(0))
    {
        return tmp<scalarField>(new scalarField(patch().size(), Zero));
    }

    const word& YName = internalField().name();

    // Initialise the input variables to the mass fractions
    scalarField psic(patchInternalField());

    const label nbrPatchi = samplePolyPatch().index();
    const fvPatch& nbrPatch = patch().boundaryMesh()[nbrPatchi];
    const fvPatchScalarField& nbrYp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(YName);
    scalarField nbrPsic(nbrYp.patchInternalField());
    mappedPatchBase::map().distribute(nbrPsic);

    switch (input_)
    {
        case none:
            FatalErrorInFunction
                << "A none input cannot be used with a non-zero transfer "
                << "coefficient" << exit(FatalError);

        case massFraction:
            // Do nothing
            break;

        case partialPressure:
            // Multiply by pressure
            {
                psic *=
                    patch().lookupPatchField<volScalarField, scalar>(pName_);

                fvPatchScalarField nbrP
                (
                    nbrPatch.lookupPatchField<volScalarField, scalar>(pName_)
                );
                mappedPatchBase::map().distribute(nbrP);
                nbrPsic *= nbrP;
            }

            // Falls through ...

        case moleFraction:
            // Convert to mole fraction
            {
                const basicSpecieMixture& mixture = composition(db());
                const scalar Wi(mixture.Wi(mixture.species()[YName]));
                const basicThermo& thermo =
                    db().lookupObject<basicThermo>(basicThermo::dictName);

                psic *= thermo.W(patch().index())/Wi;

                scalarField nbrW(thermo.W(nbrPatch.index()));
                mappedPatchBase::map().distribute(nbrW);
                nbrPsic *= nbrW/Wi;
            }
            break;
    }

    return c_*patch().magSf()*(psic - nbrPsic);
}


void Foam::semiPermeableBaffleMassFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const scalarField muEffp(turbModel.muEff(patch().index()));
    const scalarField AMuEffp(patch().magSf()*muEffp);

    valueFraction() = phip/(phip - patch().deltaCoeffs()*AMuEffp);
    refGrad() = - phiY()/AMuEffp;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::semiPermeableBaffleMassFractionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    mappedPatchBase::write(os);
    writeEntryIfDifferent<scalar>(os, "c", scalar(0), c_);
    writeEntry(os, "input", inputNames_[input_]);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        semiPermeableBaffleMassFractionFvPatchScalarField
    );
}

// ************************************************************************* //
