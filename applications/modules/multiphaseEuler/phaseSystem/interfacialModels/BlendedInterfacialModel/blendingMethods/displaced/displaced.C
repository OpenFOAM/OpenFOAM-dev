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

#include "displaced.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(displaced, 0);
    addToRunTimeSelectionTable(blendingMethod, displaced, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::displaced::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    return blendingMethod::fContinuous(blending_(), alphas, index);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::displaced::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    return constant(alphas, displacingPhasei == phasei_ ? 1 : 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::displaced::displaced
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    phasei_(interface.fluid().phases()[dict.lookup<word>("phase")].index()),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::displaced::~displaced()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::displaced::canBeContinuous(const label index) const
{
    return blending_->canBeContinuous(index);
}


bool Foam::blendingMethods::displaced::canSegregate() const
{
    return blending_->canSegregate();
}


bool Foam::blendingMethods::displaced::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return displacingPhasei == phasei_;
}


bool Foam::blendingMethods::displaced::functionOfAlphas() const
{
    return blending_->functionOfAlphas();
}


// ************************************************************************* //
