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

#include "continuous.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(continuous, 0);
    addToRunTimeSelectionTable(blendingMethod, continuous, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::continuous::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    return constant(alphas, index == index_ ? 1 : 0);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::continuous::fSegregated
(
    const UPtrList<const volScalarField>& alphas
) const
{
    return constant(alphas, 0);
}


Foam::tmp<Foam::volScalarField>
Foam::blendingMethods::continuous::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    if (!displaced_ || otherIndex_ == -1) return constant(alphas, 0);

    const dimensionedScalar& residualAlpha =
        interface_.fluid().phases()[displacingPhasei].residualAlpha();

    tmp<volScalarField> talphaSystem = alpha(alphas, -1, false);
    const volScalarField& alphaSystem = talphaSystem();

    // This is the limit of linear/hyperbolic methods as the minimum continuous
    // alpha parameters tend to zero for the continuous phase, and to one for
    // the non-continuous phases

    return
        max(alphas[displacingPhasei], residualAlpha)
       /max(1 - alphaSystem, residualAlpha)
       *pos(alpha(alphas, otherIndex_, false) - sqr(alphaSystem));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::continuous::continuous
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    phase_
    (
        interface_.fluid().phases()
        [
            dict.lookupBackwardsCompatible<word>({"continuousPhase", "phase"})
        ]
    ),
    displaced_
    (
        allowDisplaced
      ? dict.lookupOrDefault<bool>("displaced", true)
      : false
    ),
    index_(interface.contains(phase_) ? interface.index(phase_) : -1),
    otherIndex_(index_ != -1 ? interface.otherIndex(phase_) : -1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::continuous::~continuous()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::continuous::canBeContinuous(const label index) const
{
    return index == index_;
}


bool Foam::blendingMethods::continuous::canSegregate() const
{
    return false;
}


bool Foam::blendingMethods::continuous::isDisplacedBy
(
    const label displacingPhasei
) const
{
    return displaced_;
}


bool Foam::blendingMethods::continuous::functionOfAlphas() const
{
    return true;
}


// ************************************************************************* //
