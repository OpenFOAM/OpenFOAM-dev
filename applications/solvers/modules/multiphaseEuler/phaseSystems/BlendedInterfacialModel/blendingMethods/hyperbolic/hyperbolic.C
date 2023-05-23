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

#include "hyperbolic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(hyperbolic, 0);
    addToRunTimeSelectionTable(blendingMethod, hyperbolic, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label phaseSet,
    const label systemSet
) const
{
    tmp<volScalarField> x = this->x(alphas, phaseSet, systemSet);
    tmp<volScalarField> a = parameter(alphas, phaseSet, minContinuousAlpha_);
    return (1 + tanh((4/transitionAlphaScale_.value)*(x - a)))/2;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::hyperbolic
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    blendingMethod(dict, interface),
    minContinuousAlpha_
    (
        readParameters("minContinuousAlpha", dict, interface, {0, 1}, true)
    ),
    transitionAlphaScale_
    (
        readParameter("transitionAlphaScale", dict, {0, vGreat}, false)
    )
{
    if
    (
        canBeContinuous(0)
     && canBeContinuous(1)
     && minContinuousAlpha_[0].value + minContinuousAlpha_[1].value
      < 1 - rootSmall
    )
    {
        FatalErrorInFunction
            << typeName.capitalise() << " blending function for interface "
            << interface.name() << " is invalid in that it creates negative "
            << "coefficients for sub-modelled values. A valid function will "
            << "have minimum continuous alphas that sum one or greater."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::~hyperbolic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::hyperbolic::canBeContinuous(const label index) const
{
    return minContinuousAlpha_[index].valid;
}


bool Foam::blendingMethods::hyperbolic::canSegregate() const
{
    return
        canBeContinuous(0)
     && canBeContinuous(1)
     && minContinuousAlpha_[0].value + minContinuousAlpha_[1].value
      > 1 + rootSmall;
}


// ************************************************************************* //
