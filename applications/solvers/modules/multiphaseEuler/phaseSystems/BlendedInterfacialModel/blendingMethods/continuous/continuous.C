/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::continuous::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label phaseSet,
    const label systemSet
) const
{
    return
        constant
        (
            alphas,
            interface_.contains(phase_)
         && (0b01 << interface_.index(phase_) & phaseSet)
          ? 1
          : 0
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::continuous::continuous
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    blendingMethod(dict, interface),
    phase_(interface_.fluid().phases()[dict.lookup<word>("phase")])
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::continuous::~continuous()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::continuous::canBeContinuous(const label index) const
{
    return interface_.contains(phase_) && interface_.index(phase_) == index;
}


bool Foam::blendingMethods::continuous::canSegregate() const
{
    return false;
}


// ************************************************************************* //
