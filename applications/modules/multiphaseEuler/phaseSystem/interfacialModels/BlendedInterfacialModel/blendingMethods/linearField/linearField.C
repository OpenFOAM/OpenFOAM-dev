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

#include "linearField.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(linearField, 0);
    addToRunTimeSelectionTable(blendingMethod, linearField, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Pair<Foam::dimensionedScalar>&
Foam::blendingMethods::linearField::fieldValues() const
{
    if (fieldValuesPtr_.empty())
    {
        const volScalarField& field =
            interface_.mesh().lookupObject<volScalarField>(fieldName_);

        Pair<scalar> values =
            dictPtr_().lookup<Pair<scalar>>("values", field.dimensions());

        if (values.first() > values.second())
        {
            Foam::Swap(values.first(), values.second());
        }

        fieldValuesPtr_.set
        (
            new Pair<dimensionedScalar>
            (
                dimensionedScalar(field.dimensions(), values.first()),
                dimensionedScalar(field.dimensions(), values.second())
            )
        );

        dictPtr_.clear();
    }

    return fieldValuesPtr_();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linearField::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    const volScalarField& field =
        alphas[0].mesh().lookupObject<volScalarField>(fieldName_);

    const dimensionedScalar& p = fieldValues().first();
    const dimensionedScalar& f = fieldValues().second();

    const volScalarField fraction(min(max((field - p)/(f - p), zero()), one()));

    return
        (1 - fraction)*blendingMethod::fContinuous
        (
            belowBlending_(),
            alphas,
            index
        )
      + fraction*blendingMethod::fContinuous
        (
            aboveBlending_(),
            alphas,
            index
        );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linearField::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    const volScalarField& field =
        alphas[0].mesh().lookupObject<volScalarField>(fieldName_);

    const dimensionedScalar& p = fieldValues().first();
    const dimensionedScalar& f = fieldValues().second();

    const volScalarField fraction(min(max((field - p)/(f - p), zero()), one()));

    return
        (1 - fraction)*blendingMethod::fDisplaced
        (
            belowBlending_(),
            alphas,
            displacingPhasei
        )
      + fraction*blendingMethod::fDisplaced
        (
            aboveBlending_(),
            alphas,
            displacingPhasei
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linearField::linearField
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    fieldName_(dict.lookup<word>("field")),
    dictPtr_(dict.clone()),
    fieldValuesPtr_(nullptr),
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

Foam::blendingMethods::linearField::~linearField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::linearField::canBeContinuous
(
    const label index
) const
{
    return
        belowBlending_->canBeContinuous(index)
     || aboveBlending_->canBeContinuous(index);
}


bool Foam::blendingMethods::linearField::canSegregate() const
{
    return
        belowBlending_->canSegregate()
     || aboveBlending_->canSegregate();
}


bool Foam::blendingMethods::linearField::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return
        belowBlending_->isDisplacedBy(displacingPhasei)
     || aboveBlending_->isDisplacedBy(displacingPhasei);
}


bool Foam::blendingMethods::linearField::functionOfAlphas() const
{
    return
        belowBlending_->functionOfAlphas()
     && aboveBlending_->functionOfAlphas()
     && IOobject::member(fieldName_) == "alpha"
     && interface_.fluid().phases().found(IOobject::group(fieldName_));
}


// ************************************************************************* //
