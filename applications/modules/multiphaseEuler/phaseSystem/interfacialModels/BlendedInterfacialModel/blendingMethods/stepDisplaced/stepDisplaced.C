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

#include "stepDisplaced.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(stepDisplaced, 0);
    addToRunTimeSelectionTable(blendingMethod, stepDisplaced, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField>
Foam::blendingMethods::stepDisplaced::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    return blendingMethod::fContinuous(blending_(), alphas, index);
}


Foam::tmp<Foam::volScalarField>
Foam::blendingMethods::stepDisplaced::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    if (!isDisplacedBy(displacingPhasei)) return constant(alphas, 0);

    const scalar a = displacingAlphas_[displacingPhasei].value;

    // Quick optimisation for a single displacing phase
    if (displacingPhases_.size() == 1)
    {
        return pos(alphas[displacingPhasei] - a);
    }

    // Solution for multiple displacing phases
    tmp<volScalarField> integralF = posPart(alphas[displacingPhasei] - a);
    tmp<volScalarField> sumIntegralFbyF = integralF().clone();
    forAll(displacingPhases_, i)
    {
        if (displacingPhases_[i] == displacingPhasei) continue;

        const scalar a = displacingAlphas_[displacingPhases_[i]].value;

        sumIntegralFbyF.ref() += posPart(alphas[displacingPhases_[i]] - a);
    }
    return integralF/max(sumIntegralFbyF, vSmall);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::stepDisplaced::stepDisplaced
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    displacingAlphas_(interface.fluid().phases().size(), {false, NaN}),
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

        displacingAlphas_[phasei] =
            readParameter
            (
                IOobject::groupName("displacingAlpha", phase.name()),
                dict,
                {0, 1},
                NaN
            );

        if (displacingAlphas_[phasei].specified)
        {
            displacingPhases_.append(phasei);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::stepDisplaced::~stepDisplaced()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::stepDisplaced::canBeContinuous
(
    const label index
) const
{
    return blending_->canBeContinuous(index);
}


bool Foam::blendingMethods::stepDisplaced::canSegregate() const
{
    return blending_->canSegregate();
}


bool Foam::blendingMethods::stepDisplaced::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return displacingAlphas_[displacingPhasei].specified;
}


bool Foam::blendingMethods::stepDisplaced::functionOfAlphas() const
{
    return blending_->functionOfAlphas();
}


// ************************************************************************* //
