/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "stepField.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(stepField, 0);
    addToRunTimeSelectionTable(blendingMethod, stepField, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::dimensionedScalar&
Foam::blendingMethods::stepField::fieldValue() const
{
    if (fieldValuePtr_.empty())
    {
        const volScalarField& field =
            interface_.mesh().lookupObject<volScalarField>(fieldName_);

        fieldValuePtr_.set
        (
            new dimensionedScalar("value", field.dimensions(), dictPtr_())
        );

        dictPtr_.clear();
    }

    return fieldValuePtr_();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::stepField::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    const volScalarField& field =
        alphas[0].mesh().lookupObject<volScalarField>(fieldName_);

    return
        neg0(field - fieldValue())*blendingMethod::fContinuous
        (
            belowBlending_(),
            alphas,
            index
        )
      + pos(field - fieldValue())*blendingMethod::fContinuous
        (
            aboveBlending_(),
            alphas,
            index
        );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::stepField::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    const volScalarField& field =
        alphas[0].mesh().lookupObject<volScalarField>(fieldName_);

    return
        neg0(field - fieldValue())*blendingMethod::fDisplaced
        (
            belowBlending_(),
            alphas,
            displacingPhasei
        )
      + pos(field - fieldValue())*blendingMethod::fDisplaced
        (
            aboveBlending_(),
            alphas,
            displacingPhasei
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::stepField::stepField
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    fieldName_(dict.lookup<word>("field")),
    dictPtr_(dict.clone()),
    fieldValuePtr_(nullptr),
    belowBlending_
    (
        blendingMethod::New
        (
            "below",
            dict.subDict("below"),
            interface,
            allowDisplaced
        )
    ),
    aboveBlending_
    (
        blendingMethod::New
        (
            "above",
            dict.subDict("above"),
            interface,
            allowDisplaced
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::stepField::~stepField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::stepField::canBeContinuous(const label index) const
{
    return
        belowBlending_->canBeContinuous(index)
     || aboveBlending_->canBeContinuous(index);
}


bool Foam::blendingMethods::stepField::canSegregate() const
{
    return
        belowBlending_->canSegregate()
     || aboveBlending_->canSegregate();
}


bool Foam::blendingMethods::stepField::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return
        belowBlending_->isDisplacedBy(displacingPhasei)
     || aboveBlending_->isDisplacedBy(displacingPhasei);
}


bool Foam::blendingMethods::stepField::functionOfAlphas() const
{
    return
        belowBlending_->functionOfAlphas()
     && aboveBlending_->functionOfAlphas()
     && IOobject::member(fieldName_) == "alpha"
     && interface_.fluid().phases().found(IOobject::group(fieldName_));
}


// ************************************************************************* //
