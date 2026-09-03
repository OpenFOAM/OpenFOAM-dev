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

#include "linear.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(linear, 0);
    addToRunTimeSelectionTable(blendingMethod, linear, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    tmp<volScalarField> x = this->x(alphas, index);
    tmp<volScalarField> f = parameter(alphas, index, minFullyContinuousAlpha_);
    tmp<volScalarField> p = parameter(alphas, index, minPartlyContinuousAlpha_);
    return min(max((x - p())/max(f - p(), rootVSmall), zero()), one());
}


Foam::tmp<Foam::volScalarField>
Foam::blendingMethods::linear::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    if (!displaced_) return constant(alphas, 0);

    const dimensionedScalar& residualAlpha =
        interface_.fluid().phases()[displacingPhasei].residualAlpha();

    return
        max(alphas[displacingPhasei], residualAlpha)
       /max(1 - alpha(alphas, -1, false), residualAlpha)
       *(1 - fContinuous(alphas, -1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::linear
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    displaced_
    (
        allowDisplaced
      ? dict.lookupOrDefault<bool>("displaced", true)
      : false
    ),
    minFullyContinuousAlpha_
    (
        readParameters
        (
            "minFullyContinuousAlpha",
            dict,
            interface,
            {0, 1},
            1
        )
    ),
    minPartlyContinuousAlpha_
    (
        readParameters
        (
            "minPartlyContinuousAlpha",
            dict,
            interface,
            {0, 1},
            1
        )
    )
{
    forAllConstIter(phaseInterface, interface, iter)
    {
        const label i = iter.index();

        if
        (
            minFullyContinuousAlpha_[i].specified
         != minPartlyContinuousAlpha_[i].specified
        )
        {
            FatalIOErrorInFunction(dict)
                << "Both minimum fully and partly continuous alpha must be "
                << "supplied for phases that can become continuous. Only one "
                << "is supplied for " << iter().name() << "."
                << exit(FatalIOError);
        }

        if
        (
            (
                canBeContinuous(i)
             && minFullyContinuousAlpha_[i].value
             <= minPartlyContinuousAlpha_[i].value
            )
        )
        {
            FatalIOErrorInFunction(dict)
                << "The fully continuous alpha specified for " << iter().name()
                << " is not greater than the partly continuous alpha"
                << exit(FatalIOError);
        }
    }

    if
    (
        canBeContinuous(0)
     && canBeContinuous(1)
     && (
            (
                minFullyContinuousAlpha_[0].value
              + minPartlyContinuousAlpha_[1].value
              < 1 - rootSmall
            )
         || (
                minFullyContinuousAlpha_[1].value
              + minPartlyContinuousAlpha_[0].value
              < 1 - rootSmall
            )
        )
    )
    {
        FatalIOErrorInFunction(dict)
            << typeName.capitalise() << " blending function for interface "
            << interface.name() << " is invalid in that it creates negative "
            << "coefficients for sub-modelled values. A valid function will "
            << "have fully continuous alphas that are greater than one minus "
            << "the partly continuous alphas in the opposite phase."
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::~linear()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::linear::canBeContinuous(const label index) const
{
    return minFullyContinuousAlpha_[index].specified;
}


bool Foam::blendingMethods::linear::canSegregate() const
{
    return
        canBeContinuous(0)
     && canBeContinuous(1)
     && (
            (
                minFullyContinuousAlpha_[0].value
              + minPartlyContinuousAlpha_[1].value
              > 1 + rootSmall
            )
         || (
                minFullyContinuousAlpha_[1].value
              + minPartlyContinuousAlpha_[0].value
              > 1 + rootSmall
            )
        );
}


bool Foam::blendingMethods::linear::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return displaced_;
}


bool Foam::blendingMethods::linear::functionOfAlphas() const
{
    return true;
}


// ************************************************************************* //
