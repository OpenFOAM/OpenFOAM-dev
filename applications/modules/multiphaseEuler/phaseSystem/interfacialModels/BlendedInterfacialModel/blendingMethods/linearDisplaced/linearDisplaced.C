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

#include "linearDisplaced.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(linearDisplaced, 0);
    addToRunTimeSelectionTable(blendingMethod, linearDisplaced, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField>
Foam::blendingMethods::linearDisplaced::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    return blendingMethod::fContinuous(blending_(), alphas, index);
}


Foam::tmp<Foam::volScalarField>
Foam::blendingMethods::linearDisplaced::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    if (!isDisplacedBy(displacingPhasei)) return constant(alphas, 0);

    const volScalarField& x = alphas[displacingPhasei];
    const scalar f = minFullyDisplacingAlphas_[displacingPhasei].value;
    const scalar p = minPartlyDisplacingAlphas_[displacingPhasei].value;

    // Quick optimisation for just one displacing phase
    if (displacingPhases_.size() == 1)
    {
        return min(max((x - p)/(f - p), zero()), one());
    }

    // Solution for multiple displacing phases
    tmp<volScalarField> integralF =
        pos(x - p)*neg0(x - f)*sqr(x - p)/(f - p)/2
      + pos(x - f)*(x - (p + f)/2);
    tmp<volScalarField> sumIntegralFbyF;
    forAll(displacingPhases_, i)
    {
        const volScalarField& x = alphas[displacingPhases_[i]];
        const scalar f = minFullyDisplacingAlphas_[displacingPhases_[i]].value;
        const scalar p = minPartlyDisplacingAlphas_[displacingPhases_[i]].value;
        auto integralFByF =
            pos(x - p)*neg0(x - f)*(x - p)/2
          + pos(x - f)*(x - (p + f)/2);
        if (!sumIntegralFbyF.valid())
        {
            sumIntegralFbyF = eval(integralFByF);
        }
        else
        {
            sumIntegralFbyF.ref() += integralFByF;
        }
    }
    return integralF/max(sumIntegralFbyF, vSmall);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linearDisplaced::linearDisplaced
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    minFullyDisplacingAlphas_(interface.fluid().phases().size(), {false, NaN}),
    minPartlyDisplacingAlphas_(interface.fluid().phases().size(), {false, NaN}),
    displacingPhases_(),
    blending_
    (
        blendingMethod::New
        (
            "blending",
            dict.subDict("blending"),
            interface,
            false
        )
    )
{
    if (!allowDisplaced)
    {
        FatalIOErrorInFunction(dict)
            << "Blending method " << typeName << " selected as a sub-blending "
            << "in a context in which displaced blending is not allowed"
            << exit(FatalError);
    }

    forAll(interface.fluid().phases(), phasei)
    {
        const phaseModel& phase = interface.fluid().phases()[phasei];

        if (interface.contains(phase)) continue;

        minFullyDisplacingAlphas_[phasei] =
            readParameter
            (
                IOobject::groupName("minFullyDisplacingAlpha", phase.name()),
                dict,
                {0, 1},
                NaN
            );
        minPartlyDisplacingAlphas_[phasei] =
            readParameter
            (
                IOobject::groupName("minPartlyDisplacingAlpha", phase.name()),
                dict,
                {0, 1},
                NaN
            );

        if
        (
            minFullyDisplacingAlphas_[phasei].specified
         != minPartlyDisplacingAlphas_[phasei].specified
        )
        {
            FatalIOErrorInFunction(dict)
                << "Both minimum fully and partly displacing alpha must be "
                << "supplied for a displacing phases. Only one is supplied for "
                << phase.name() << "." << exit(FatalIOError);
        }

        if
        (
            minFullyDisplacingAlphas_[phasei].specified
         && minFullyDisplacingAlphas_[phasei].value
         <= minPartlyDisplacingAlphas_[phasei].value
        )
        {
            FatalIOErrorInFunction(dict)
                << "The fully displacing alpha specified for " << phase.name()
                << " is not greater than the partly continuous alpha"
                << exit(FatalIOError);
        }

        if (minFullyDisplacingAlphas_[phasei].specified)
        {
            displacingPhases_.append(phasei);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::linearDisplaced::~linearDisplaced()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::linearDisplaced::canBeContinuous
(
    const label index
) const
{
    return blending_->canBeContinuous(index);
}


bool Foam::blendingMethods::linearDisplaced::canSegregate() const
{
    return blending_->canSegregate();
}


bool Foam::blendingMethods::linearDisplaced::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return minFullyDisplacingAlphas_[displacingPhasei].specified;
}


bool Foam::blendingMethods::linearDisplaced::functionOfAlphas() const
{
    return blending_->functionOfAlphas();
}


// ************************************************************************* //
