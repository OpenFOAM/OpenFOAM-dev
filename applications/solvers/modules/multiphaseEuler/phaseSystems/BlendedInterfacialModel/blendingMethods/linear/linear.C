/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
#include "one.H"
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label phaseSet,
    const label systemSet
) const
{
    tmp<volScalarField> x = this->x(alphas, phaseSet, systemSet);
    tmp<volScalarField> f =
        parameter(alphas, phaseSet, minFullyContinuousAlpha_);
    tmp<volScalarField> p =
        parameter(alphas, phaseSet, minPartlyContinuousAlpha_);
    return min(max((x - p())/max(f - p(), rootVSmall), zero()), one());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::linear
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    blendingMethod(dict, interface),
    minFullyContinuousAlpha_
    (
        readParameters
        (
            "minFullyContinuousAlpha",
            dict,
            interface,
            {0, 1},
            true
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
            true
        )
    )
{
    forAllConstIter(phaseInterface, interface, iter)
    {
        const label i = iter.index();

        if
        (
            minFullyContinuousAlpha_[i].valid
         != minPartlyContinuousAlpha_[i].valid
        )
        {
            FatalErrorInFunction
                << "Both minimum fully and partly continuous alpha must be "
                << "supplied for phases that can become continuous. Only one "
                << "is supplied for " << iter().name() << exit(FatalError);
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
            FatalErrorInFunction
                << "The fully continuous alpha specified for " << iter().name()
                << " is not greater than the partly continuous alpha"
                << exit(FatalError);
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
        FatalErrorInFunction
            << typeName.capitalise() << " blending function for interface "
            << interface.name() << " is invalid in that it creates negative "
            << "coefficients for sub-modelled values. A valid function will "
            << "have fully continuous alphas that are greater than one minus "
            << "the partly continuous alphas in the opposite phase."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::~linear()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::linear::canBeContinuous(const label index) const
{
    return minFullyContinuousAlpha_[index].valid;
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


// ************************************************************************* //
